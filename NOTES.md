# Notes

Working notes on the notation system. The code is the reference for what is
implemented; this is for the reasoning behind it and for what is still open.

## The main idea: a notation is a temperament

A notation system is a linear map from the interval vectors of a just
intonation subgroup to **notation coordinates**

```
(octave, fifth, accidental_1, accidental_2, ...)
```

exactly like a temperament mapping. The first two coordinates are always the
octave `2/1` and the fifth `3/2`; together they give the nominals and the
sharps and flats, a sharp being seven fifths less four octaves, `2187/2048`.
Each further coordinate counts one accidental, which raises or lowers by a
small interval that is not a sharp.

Its kernel is the set of intervals it spells identically.

The real temperament always sits inside the notational one:

```
ker(notation) is contained in ker(temperament)
```

If it were not, two intervals that sound different would be written the same
way. Everything else follows from that:

- `rank(notation) >= rank(temperament)` always.
- `rank(notation) - 2` is the number of accidentals.
- `rank(notation) == rank(temperament)` makes spelling a **bijection**: every
  pitch has exactly one spelling and there is nothing to choose.
- `rank(notation) > rank(temperament)` leaves a kernel of **enharmonics**,
  spellings that name the same pitch. These are not a defect - they are what
  preserves the harmonic intent the temperament discarded. `C#` and `Db` are
  distinct in meantone and the same in 12et; 22et writing its major third as
  `D#` and 41et writing its third as `Fb` are defining properties of those
  temperaments, not mistakes.

A notation is therefore never derived automatically in the sense of "the right
answer". Each (notation, temperament) pair is a design decision, and the
library's job is to offer the choices that make sense and let the caller pick.

## How a notation is built

The data is a subgroup, a list of accidentals, and a list of notational commas,
with `accidentals + commas == dim - 2`. The mapping is derived: stack the
generators (octave, fifth, accidentals) on top of the commas into a square
matrix and invert it. `mapping * transpose(square) = [I | 0]` is the condition
that each generator gets its unit vector and each comma goes to zero, so

```
mapping = first `rank` rows of transpose(inverse(square))
```

`integer_inverse` fails unless the determinant is ±1, which is exactly the
check that the generators and commas together span the subgroup. If they do
not, some interval has no integer spelling. `Notation::assemble` is this.

## Accidentals

Derived, not given. For each prime beyond 3, shift it along the fifth chain by
0, 1, -1, 2, -2, ... fifths, octave-reduce each candidate to within a tritone
of a unison, and take the first that lands within **half an apotome** (56.84
cents, `600 * (7 * log2(3) - 11)`).

That bound is the load-bearing part: anything wider is closer to the
neighbouring point of the fifth chain, so it belongs there as a sharp or flat
instead. Bounding accidentals here is what lets every prime be written
directly, with no augmented or diminished interval needed.

The rule generalises further than the cases it was designed for:

```
5  -> 81/80      7  -> 64/63      11 -> 33/32
13 -> 1053/1024  17 -> 4131/4096  19 -> 513/512
```

The 13, 17 and 19 choices are the Helmholtz-Ellis ones, for free. `33/32` at
53.3 cents is the widest accidental in normal use, 3.6 cents inside the bound;
`27/26` at 65.3 cents is what the bound excludes for 13.

By construction an accidental has exponent ±1 on its own prime and support
`{2, 3, p}`. Normalising them all to ascending gives the conventional
directions for free: 5 and 7 land on coordinate -1 (lowered), 11 on +1.

## What decides whether a notation exists

A rank `r` notation exists exactly when some subset `S` of the primes beyond 3
with `|S| = r - 2` makes

```
[octave, fifth, accidentals of S] together with a basis of ker(temperament)
```

span the whole subgroup lattice. Equivalently, the temperament's images of
those generators generate `Z^rank(temperament)`. (Testable with
`solve_diophantine(transpose(rows), eye(dim)).is_ok()`.)

- For `r == rank(temperament)` the matrix is square and this is `det = ±1`.
- The **leftmost `rank x rank` block of the HNF mapping** is the instance of
  this test with `S =` the first `r - 2` primes. Its determinant is exactly 1
  when a rank `r` notation exists with that choice of `S`. For rank 2 there is
  no choice, so it is exact. From rank 3 up it is only sufficient: marvel keeps
  `81/80` with determinant 1 but `64/63` with determinant 2, so the leftmost
  block could say 2 while a notation exists via another subset. Cheap enough to
  search.
- `r = dim` (the just intonation notation) always spans, so **some** notation
  always exists.

## Uniqueness: where all the difficulty lives

This is the thing worth remembering.

- At `r == rank(temperament)` the notation is **forced**. `ker(N) = ker(T)`,
  the matrix inversion produces it, there is nothing to choose and no heuristic
  is needed.
- At `r > rank(temperament)` the notation is **underdetermined**. Writing
  `T = t . N`, the freedom is a map into `ker(t)` - and `ker(t)` is precisely
  the enharmonic lattice. So choosing a notation above the temperament's rank
  *is* the simplification problem. Any rule for it is a heuristic.

The enharmonic lattice is worth computing directly:

```
images[i] = temperament.map(generator[i])       // rank x r
enharmonics = kernel_left(images)               // rank r - rank(temperament)
```

An equal temperament is rank 1, so *every* notation of one is above its rank.
That is why they need special handling, and why they came first.

## Every rank, one mechanism (implemented)

`Notation::options` now works the same way at every rank, and the equal
temperament special case is gone. Each accidental is put in one of four classes:

- **tempered out** - the temperament maps it to zero, so it would raise by
  nothing. Never offered, and its comma is just itself.
- **passed over** - worth exactly what an accidental already kept is worth, up to
  direction, so it would be that accidental over again. Never offered either.
  `81/80` and `64/63` being one interval is a property of 41et, not a second
  symbol to read; likewise `33/32` and `64/63` in 31et, and `81/80` and `64/63`
  in pele, so this is not only an equal temperament rule.
- **necessary** - part of the smallest subset of what is left that makes a
  notation possible at all. Always kept.
- **optional** - the rest. The run takes them on one at a time, from keeping only
  the necessary ones to keeping them all.

**Necessary** is a search, `Notation::necessary`. A subset works when the images
of the octave, the fifth and those accidentals generate the whole tempered
lattice - otherwise some interval the temperament distinguishes has no spelling.
Accidentals offered first are kept first where there is a choice. The smallest
subset has `rank - 2` members exactly when a notation of the temperament's own
rank exists, so the old `minimal` is this search and needs no separate code.

The rank 1 tests of the old code turn out to be this search written out: the
fifth chain reaching every note is `gcd(divisions, fifth) == 1`, i.e. the empty
subset spanning. The one thing the search does not say by itself is the **ups and
downs rule: for an equal temperament the finest accidental is one step**.
Spanning only needs `gcd(divisions, fifth, steps) == 1`, which 25et meets with
`81/80` worth two steps - but a notation with no symbol for a single step can
only reach single steps by walking the fifth chain, which is not how anyone
writes such a temperament. So `Notation::candidates` gives an equal temperament
nothing to work with unless some accidental is worth one step, and puts that one
first. **25et, 51et and 54et therefore have no notation** over `2.3.5`, and the
answer for them is a subgroup whose accidental does fit: 24et is contorted over
`2.3.5` and saturates to 12et, but over `2.3.5.11` its quartertone is `33/32`,
worth exactly one step, and it notates fine. That entry is in the list now.

This is the one place where rank 1 is still singled out, and it has to be: for
rank above 1 no accidental can generate the tempered lattice by itself, so there
is nothing for "worth one step" to generalise to.

## Replacement: the one remaining choice

Each accidental that is not kept gets **one comma, fixed once**, and the kernel
of a notation is spanned by the commas of the accidentals it dropped. Taking an
accidental on only ever removes a comma, so **kernels nest** by construction - a
smaller notation spells alike everything a larger one does, and a larger
notation's spelling can always be simplified onto a smaller one's. This is
checked over the whole list and every equal temperament to 72 by
`cargo run --example verify`.

It is also what makes unimodularity automatic for every subset. A comma may be
built from the octave, the fifth, the necessary accidentals, and the optional
ones the run has already taken on - exactly the ones kept wherever this
accidental is dropped. One that is never kept may use every optional accidental.
Nothing is ever built from an accidental that is never kept, so substituting the
commas back is triangular and terminates at the fifth chain. That last clause
matters: letting a never-kept accidental stand in for another one only stacks one
substitution on top of another, and it is why 41et's `7/4` is `vBb` - `64/63`
reaching for `81/80`, which the notation does keep - rather than a detour on the
chain.

The comma is `accidental - replacement`, and the replacement has to be worth
what the accidental is worth. That usually leaves a choice, settled in two
steps.

1. **A plain stack of the accidentals still available**, if there is one, with
   no octaves and no fifths. The prime then keeps the nominal and the sharps
   that just intonation gives it and only the accidental count changes. This is
   what writes `33/32` as two syntonic commas in 41et and as one septimal comma
   in 31et, and it reproduces every spelling the old equal temperament rule
   produced. Where several stacks work, the lower primes are preferred, which
   is the only tie-break needed.
2. **Otherwise the simplest comma.** The fifth chain has to be walked and how
   far is a choice again; the replacements that work differ by the commas the
   notation could temper out, so they form a coset of that lattice, and the
   element of smallest `Weighting::Wilson` norm is taken (LLL then `cvp_exact`).

Step 2 is the part that is not forced, and the case for it is thinner than it
first looked. The alternative is to extend step 1's preference order to the
octave and the fifth - avoid the fifth first, then the octave, then the higher
accidentals. That runs away, because avoiding the fifth at any price is what
lexicographic preference means: asked to replace a one step accidental where only
a two step one is available, it answers with a stack of **25** of them and an
octave off rather than take five fifths. But the ups and downs rule heads that
off at rank 1 by never keeping a coarse accidental while a finer one is dropped,
and with both rules in place the two tier 2 choices agree on **every** equal
temperament to 99 and everywhere on the temperament list but one: 7-limit magic,
where the simplest comma gives `225/224` and `vvA#` and the preference order
gives `864/875` and `^^^Bbb`. Two marks and a comma the temperament is named for,
against three marks and a comma nobody would name.

So step 2 is currently worth one better answer and some insurance. The insurance
is the point: nothing at rank 2 and up plays the part the ups and downs rule
plays at rank 1, so nothing there stops the preference order reaching for a long
stack, and the list is thirty entries. Dropping step 2 would make the whole
derivation metric-free, which is worth wanting - if a wider search of rank 2
temperaments never makes the two disagree by more than magic does, drop it.

Using the simplest comma *everywhere*, for the record, is clearly wrong: it
spells 41et's `7/4` as `vvA#` instead of `vBb` and 31et's `11/8` as `vvGb`
instead of `^F`, because a simpler comma is not the same thing as the intended
one. The stack has to come first.

Wilson weighting is a metric rather than a tuned weight, and
`reduced_comma_basis` already uses it, so this is not the weighted cost function
the dead end above was about. It also cannot break kernel nesting the way that
one did, since nesting no longer depends on how the comma is chosen.

### What the run looks like now

Every equal temperament comes out exactly as the old code had it, which is the
main evidence that the general rules are the right ones. What is new is the
middle of the run for rank 2 and up.

- **huygens** gains a rank 3: keep `64/63`, write `11/8` as `vvGb`.
- **11-limit marvel and miracle** gain a rank 4: keep `81/80` and `64/63`, write
  `11/8` as `^>F`, one of each, since `33/32` is worth exactly what the two are
  worth together.
- **7-limit miracle, orwell and magic** gain a rank 3, keeping only `81/80` and
  writing `7/4` as `vvA#` over the comma `225/224`, which all three temper out.
  The old code refused to guess here and jumped to the next notation up.
- **pele loses** its rank 4. It tempers out `5120/5103`, so its accidentals for 5
  and 7 are one interval and a notation with both is the same notation twice.

Spellings in the bottom notation of a large equal temperament are wild - 41et in
the 11-limit writes `11/8` as `D###` where the old code gave `Abbb` - because
the substitution compounds: `33/32` becomes two syntonic commas and each of
those becomes twelve fifths. That is the price of kernel nesting, and it is
`simplify`'s business, not the mapping's.

## Loose ends and known limits

- **`simplify` is not started.** The pieces it needs are in place: the
  enharmonic lattice is `kernel_left` of the generator images, and reduction is
  a small closest-vector problem in a metric that ignores the octave
  coordinate. `diophantine` has `lll`, `nearest_plane` and `cvp_exact`; a
  bounded box search after reduction also works and takes an arbitrary cost
  function, which matters because the natural cost is not quadratic (an
  asymmetric nominal window, accidentals weighted against sharps). It now has
  real work waiting for it: the bottom notation of a large equal temperament
  spells `11/8` as `D###`, and kernel nesting guarantees such a spelling can be
  reduced onto the one the notation above it gives.
- **Derived `PartialEq` on `Notation` compares the stored comma basis**, so two
  notations with the same mapping built from different bases of the same kernel
  compare unequal. Nothing relies on it any more now that `options` builds the
  whole run from one set of commas, but it is still a trap.
- **Accidental symbols are positional**, so the same comma prints as `^` in one
  notation of a temperament and `>` in another. Fine for debugging output;
  wrong if symbols should be stable across a temperament's run, in which case
  key them to the prime.
- **Four symbols means 13-limit is the ceiling** for `from_ji`. Nothing else
  limits the subgroup.
- **`Temperament::et` silently defactors.** `et(24, 2.3.5)` saturates to
  `[12, 19, 28]`, i.e. to 12et, because the 24et patent mapping of the 5-limit
  is contorted. 24et as such cannot be asked for this way.
- **An accidental defined as one step** - what ups-and-downs uses in general -
  is still not available, so an equal temperament with no `(2, 3, p)` accidental
  worth one step gets no notation rather than an arbitrary one. Over `2.3.5`
  that is 25, 51 and 54 up to 72; over `2.3.5.7`, 25, 54 and 57. A wider
  subgroup is the answer for these, and the `24et` entry shows it working.
- **`Weighting` is not exposed on `options`.** The fallback in step 2 hardcodes
  the default. Tenney and Wilson agree exactly - same commas, same spellings -
  over the temperament list and over every equal temperament to 99 in the 5, 7
  and 11-limit, so nothing turns on it yet, and the fallback may well go away
  entirely.
- **Whether `225/224` may be an accidental** was left open. Nothing enforces
  the `{2, 3, p}` support rule, because accidentals are derived rather than
  given, and the derivation cannot produce anything else.
