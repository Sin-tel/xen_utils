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

**Necessary** is a search, `Search::necessary`. A subset works when the images
of the octave, the fifth and those accidentals generate the whole tempered
lattice - otherwise some interval the temperament distinguishes has no spelling.
Accidentals offered first are kept first where there is a choice. The smallest
subset has `rank - 2` members exactly when a notation of the temperament's own
rank exists, so the old `minimal` end of the run is this search and needs no
separate code.

The rank 1 tests of the old code turn out to be this search written out: the
fifth chain reaching every note is `gcd(divisions, fifth) == 1`, i.e. the empty
subset spanning. The one thing the search does not say by itself is the **ups and
downs rule: for an equal temperament the finest accidental is one step**.
Spanning only needs `gcd(divisions, fifth, steps) == 1`, which 25et meets with
`81/80` worth two steps - but a notation with no symbol for a single step can
only reach single steps by walking the fifth chain, which is not how anyone
writes such a temperament. So `Search::candidates` gives an equal temperament
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
   produced. Where several stacks work, **the shortest one** is taken - which
   in the notation basis is the comma worth the fewest marks. See the notation
   basis section below: this step is exactly "stay inside the accidentals", and
   its tie-break is exactly "write fewer symbols".
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
derivation metric-free, which is worth wanting - **two metric-free replacements
for it were tried in the notation basis and both are worse.** See below.

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

## The notation basis

The just intonation notation - octave, fifth, one accidental per prime beyond 3
- is a **change of basis and nothing more**. Its generator matrix is triangular
with `+-1` down the diagonal, since each accidental has exponent `+-1` on its
own prime and support `{2, 3, p}`, so it is unimodular and `to_notation` and
`to_just` are mutually inverse there. Any comma may be read in either basis
without losing anything. `cargo run --example basis` prints both and asserts the
round trip.

What it buys is legibility. `5120/5103` is `[10, -6, 1, -1]` over the primes and
`[0, 0, -1, 1]` in the notation basis, which says `81/80` and `64/63` are one
interval and says it on sight.

**Counting marks needs no matrix.** No generator other than the accidental for
`p` touches `p` at all, so the notation coordinate of that accidental is `+-` the
exponent of `p`. The number of accidental marks an interval is worth is
therefore `sum |e_i|` over the primes beyond 3, read straight off the prime
coordinates.

### What the basis makes obvious

Put the comma lattice in notation coordinates, columns `[octave, fifth, a_1,
...]`, and take the Hermite normal form. It splits by **where the pivot lands**:

- pivot in an accidental column, and the row has no octave and no fifth to it -
  a notational comma of the step 1 kind;
- pivot in the octave or fifth column - an **enharmonic**.

That is math.md's exact sequence `0 -> ker N -> ker T -> E -> 0` made
computational and canonical. 22et comes out as

```
[0, 0, 0, 1]      64/63 is tempered out
[1, 0, -22, 0]    the octave is 22 ups
[0, 1, -13, 0]    the fifth is 13 ups
```

which is the whole of ups and downs in 22et, read off a normal form. Every comma
the search picks that has no octave and no fifth to it agrees with the matching
row up to sign, over the whole list.

The three classes `Search` sorts accidentals into are the same rows seen from
the other side. Writing `K_0` for the sublattice with octave and fifth both
zero - the relations among accidentals alone - an accidental is **tempered out**
when `e_i` is in `K_0`, **passed over** when `e_i +- e_j` is, and gets a step 1
stack when some other element of `K_0` has `+-1` on it.

### What it does not fix

Step 2 is where the fifth chain has to be walked, and the basis does not settle
it. Two metric-free rules were implemented and swept over the temperament list
and every equal temperament to 99 across six subgroups (912 notations,
`cargo run --example dump`, which exists to be diffed):

- **Fewest accidental marks.** Runs away along the fifth chain: huygens writes
  `11/8` as `E##`, eighteen fifths up, to save two marks. 37 lines move, and
  19et's recommendation falls from `vBb` to `Bbb`.
- **Fewest symbols**, an unweighted `L1` on the notation coordinates bar the
  octave, so one fifth of walk costs one mark. 64 lines move: 12et's `11/8` goes
  `F#` to `Gb`, and 9et, 11et and 18et all shift.

The reason is structural. Marks and fifth-walk are two competing costs; any
lexicographic order on them runs away in one direction (the preference order
runs away into accidentals, fewest-marks runs away into fifths), and an
unweighted sum is an empirical claim that one fifth costs one mark. `sopfr`
balances them because it is a single norm on the actual interval, which is
exactly what math.md section 6 says the notation coordinates cannot carry: there
is no prime on the fifth axis. **Step 2 stays.**

### What it did fix

Step 1's tie-break. It used to be a lexicographic preference over the available
accidentals, lower primes first, with a nearest-to-zero rule inside the
arithmetic progression where the count was not pinned - a recursive function
with three separate conventions in it. Replacing the whole of it with **the
shortest stack** loses nothing and gains something:

- Dropping the preference order entirely, for any solution at all, changes
  **17 of 912** lines. Every one is the same shape: the `33/32` comma of an
  11-limit notation keeping both `81/80` and `64/63`. Everywhere else
  availability already pins the stack up to sign, and the old machinery was
  deciding nothing.
- Taking the shortest stack instead changes **16** lines, the same shape again,
  and fifteen of them strictly reduce the mark count. 72et is the one to look
  at: `11/8` was `^^^F`, three syntonic commas, and is now `^>F`, one syntonic
  and one septimal, since 72et makes `33/32` worth one plus two steps as well as
  three. That is the reading 11-limit marvel and miracle already had.

So the one case where step 1 ever had a choice is now settled by counting
symbols rather than by preferring low primes, which is both shorter to state and
better where it differs. `a_stack_is_the_shortest_one_that_works` pins 72et
down.

### Sweeping equal temperaments proves less than it looks

An equal temperament sweep is a **coverage test, not a judgement**. Temperaments
can be arbitrarily bad and most of the ones a sweep turns up are nobody's
notation, so the sweep is worth having for panics and for diffing a rule change
and worth nothing for deciding whether an answer is sensible. The interesting
cases are rank 2 and up.

Counting the same fifth chain walks with rank 1 and the higher ranks apart makes
the difference plain. Over `data/temperaments_big.txt` plus every equal
temperament to 99, of 43 commas from rank 2 and up:

- **25 never touch the chain at all**, so they are step 1;
- the other **18 land on 5, 7, 12 or 19 fifths** - the limma, the apotome, the
  pythagorean comma and the next convergent - and nowhere else.

The odd rows in the table, 2 fifths at a whole tone and 3 at a minor third and
26 at a major second, are **equal temperaments only**. So the intuition that a
step 2 comma always respells a prime onto a near-closure of the fifth chain is
exactly right where it matters, and the noise was coming from the sweep.

That column is thin - 43 commas - and thin is the honest state of the evidence.
Lengthening the named list, especially at rank 2 and rank 3, is what would fill
it in.

## Which notation to recommend

`from_temperament` no longer takes an end of the run. It returns **the smallest
notation that keeps every nominal**: one where each prime is spelled on the
letter just intonation gives it. Each further accidental is another symbol to
read, so smaller is better - but not at the price of moving a prime onto another
letter.

41et is the case that makes it obvious once the run is laid out:

```
   [2]  5/4 Fb5   7/4 Cbb6   11/8 D###5
-> [3]  5/4 vE5   7/4 vBb5   11/8 ^^F5
   [4]  5/4 vE5   7/4 vBb5   11/8 >F5
```

`[2]` has all three primes off their nominals; `[4]` only turns the two marks of
`[3]` into one of another kind, so `[3]` is the one wanted.

**Sharps and flats do not count.** Seven fifths leave the letter alone, so the
test is on the fifth coordinate mod 7, not on the coordinate itself. Flattone is
the case that forces this:

```
-> [2]  5/4 E5  11/8 F#5
   [3]  5/4 E5  11/8 ^F5
```

`F#` is not `^F`, but it is still an `F`, so `[2]` qualifies and there is no
reason to take on an accidental for 11.

The nominal just intonation gives a prime needs no separate notation to compute.
The accidental `a` has exponent `s = +-1` on the prime, so the prime sits
`-s * a[1]` fifths along the chain - `just_nominal` is that one line. The
syntonic comma has `s = -1` and four threes, hence 5 on `E`.

Which tier of `comma` was used is what decides this in practice: a replacement
built as a plain stack of accidentals keeps the nominal by construction, and only
a walk along the fifth chain can move it. Flattone shows a walk landing on the
same letter anyway, so the test has to be on the spelling rather than on which
tier ran.

**Where nothing qualifies** the recommendation falls back to the last notation.
That happens when every accidental has to be replaced by a chain walk that misses
- 13et over `2.3.5` is the smallest, its syntonic comma being worth two steps
while its fifth chain still reaches every note, so it escapes the refusal that
25et gets. Over every equal temperament to 200 across fourteen subgroups, 1157 of
2245 runs have no notation keeping its nominals, and **not one of them offers a
choice**: every single one is a run of length one. So the fallback never actually
decides anything, and the rule is unambiguous wherever there is anything to
decide. `a_run_with_a_choice_always_keeps_its_nominals_somewhere` pins that down.

## Enharmonics

`ker(notation)` is what the notation spells alike; the **enharmonic lattice** is
what the temperament calls a unison and the notation still spells apart. Between
them they account for everything the temperament tempers out, which `verify` now
checks:

```
rank(ker notation) + rank(enharmonics) == dim - rank(temperament)
```

`Notation::enharmonics(temperament)` is `kernel_left` of the generator images,
returned as interval vectors so that it matches `commas()`. `to_notation` puts
one back in notation coordinates, and that is the form worth reading: 22et's
rank 3 notation gives `[0, 1, -13]` and `[1, 0, -22]`, i.e. its fifth is thirteen
ups and its octave twenty two - the whole of ups and downs in 22et, which the
mapping always implied but never stated.

Every notation of an equal temperament has one, and has to: `assemble` sends each
generator to its own unit vector, so a notation is free on the octave, the fifth
and its accidentals and can never close the circle of fifths. 12et's enharmonic
is the pythagorean comma, 7et's the apotome, 5et's the limma. That is the
`C# != Db` property, not a defect.

## Where the code lives

- `notation.rs` - the `Notation` type: its mapping, its kernel, its enharmonics,
  `assemble`, `note`, and the derivation of a single accidental.
- `search.rs` - everything behind `Notation::options`. `Search` holds the
  temperament, the accidentals and their images; `Plan` is how it sorts them into
  necessary, optional and never kept. `shortest_stack` is step 1 of the comma
  choice and `Search::simplest_comma` is step 2. Tested through `options`.
- `simplify.rs` - `Simplifier`: the comma lattice reduced once, then a seeded
  walk per interval. `cargo run --example simplify` walks 41et.
- `cargo run --example trace` - the derivation of one temperament, a step at a
  time, with every replacement that works listed and the one taken marked.
  Defaults to 72et and takes any name in `data/temperaments_big.txt`. It reads
  the answer back out of `Notation::options` rather than re-deriving it, and
  asserts that everything it lists is worth what the accidental is worth, so it
  cannot quietly drift from the code.
- `util.rs` - integer vector helpers, with no music in them.

## Simplifying

Which interval a tempered pitch "is" is a closest vector problem on the
temperament's comma lattice, and three things about it took working out.

**Which lattice.** Not the enharmonic lattice - that stays inside the lattice the
notation's generators span, which for 41et's recommended notation is only
`{2, 3, 5}`, so it can never answer `7/4` and says `225/128` instead. The whole
comma lattice is what is wanted, and `commas + enharmonics = ker(temperament)`
exactly, checked over 636 notations. So simplifying depends on the temperament
alone and the notation only spells the answer.

**Which basis.** Wilson and Tenney weight *prime exponents*, so the search has to
happen in just intonation coordinates. Weighting notation coordinates instead is
meaningless: there is no prime on the fifth axis. That is also why `to_notation`
and `to_just` are named as they are now - both bases are in play at once and
"interval" for the notation side had become unreadable.

**Which norm.** Wilson is an **L1** norm - `sopfr`, the sum of the prime factors
of numerator and denominator. The quadratic form the lattice algorithms take is
only a proxy for it, and the two disagree: under L2 a step of 41et is `45/44`,
spelled `vvvC#`, where under L1 it is `64/63` at `^C`. So LLL and `cvp_exact`
seed the search and a bounded walk finishes it under `sopfr`. Ties are settled by
the quadratic norm and then by the coordinates, so the answer never depends on
the order the ball is walked in.

What comes out for 41et is the interval table anyone would write by hand: `8/7`,
`7/6`, `6/5`, `5/4`, `9/7`, `4/3`, `7/5`, `3/2`, `8/5`, `5/3`, `7/4`, `15/8`.

**The answer is the first of a list.** The walk passes over every interval in
every ball it steps through, so keeping them costs nothing and `candidates`
returns them sorted - 625 to 750 of them for an equal temperament of the 11
limit. That is a neighbourhood rather than the whole coset, so only the right end
of the list means anything, but the right end is what an end user cycles through.
`cargo run --example candidates` shows the four steps of 41et where the simplest
reading is not the most convenient spelling.

One thing to decide there: candidates are distinct *intervals*, and several of
them can share a spelling - at 14 steps both `81/64` and `80/63` are `E5`, since
they differ by a comma the notation tempers out. Right for cycling through
readings, noise for cycling through spellings.

Two things were tried and are settled.

**Leaving the octave out of the cost is wrong.** Free powers of two let the
search buy any number of them to avoid one high prime, which inverts exactly the
identifications wanted: 41et's `14/11` becomes `81/64`, `9/7` becomes `32/25`,
`7/5` becomes `45/32` and `11/7` becomes `128/81`. Octave equivalence is a
reasonable thing to want elsewhere, but `sopfr` is a statement about the ratio
and the 2s are part of the ratio.

**`nearest_plane` would do as a seed, and would not help.** Its seed differs from
`cvp_exact`'s on 4700 of 59109 targets, and the walk converges to the same answer
every time for anything within nineteen generators of the unison - the ten that
differ are the same twenty-octaves-out stalls the radius has. But the seed is
1.6% of the running time and the walk is the rest, so swapping it saves 0.6% of
nothing. Worth remembering if the lattice ever gets big enough for `cvp_exact`'s
enumeration to bite, since the walk does not need an exact seed.

## Performance

This is meant for a DAW's UI thread, not its audio thread, so constructing a
`Notation` or a `Simplifier` is not a concern - built once, reused after.
`cargo run --release --example bench` times what does run repeatedly:
`to_notation`, `to_just`, `simplify` and `candidates`.

`to_notation` and `to_just` are 35-40 nanoseconds regardless of subgroup - one
allocation, one pass over a handful of coordinates. Not worth a second thought.

`simplify` and `candidates` are the real cost, since `Simplifier::search` walks
a ball of `width.pow(lattice rank)` points, each a candidate interval - 625 of
them for 41et's 11-limit notation, 3125 for the 13-limit. It used to compute
the sort key up to three times per point, each one cloning the whole
coordinate vector just to break ties with it - the tie-break the comment
already promised is exactly what sorting `(rank, candidate)` pairs gives for
free, no separate copy needed. Fixed by ranking each candidate once, with a
key that carries no coordinates, and building each candidate as a single
`clone` of `best` mutated in place instead of a `combination` and a `subtract`
composed from scratch. Six allocations a point down to one: 41et's 11-limit
notation went from 165 to 55 microseconds a call, the 13-limit from 820 to 257.

That remaining allocation is one clone of the interval per point examined, and
is not coming out: every point may end up in the returned list, so it has to
own its coordinates. Nothing about `Matrix`'s nested-`Vec` representation is
the bottleneck here - `to_notation` and `to_just` already show that reading a
small matrix is noise next to a single allocation, and `search`'s cost is the
combinatorial walk, not how the lattice it reads is laid out. A flatter,
strided matrix type is a reasonable thing to want for other reasons, but it
would not move these numbers.

At a few tens to a few hundred microseconds a call, simplifying every visible
note every redraw is affordable at any subgroup this library is likely to see
in practice. If a future subgroup ever makes the lattice rank large enough for
this to matter again, the width of the ball - not its allocation pattern - is
where the exponential lives.

## Loose ends and known limits

- **The simplifier is a local search, not a proof.** `Simplifier` reduces the
  temperament's comma lattice once, takes the `cvp_exact` answer as a seed, and
  then walks it under the norm actually wanted. It agrees with a wider step for
  every interval within nineteen generators of the unison, over 59109 swept; the
  stalls past that are twenty octaves out, where the norm is minimized by
  trading a pile of one prime for a pile of another. See `SEARCH_RADIUS`.
- **`Notation` has no `PartialEq`.** It used to derive one, comparing the
  stored comma basis, which meant two notations with the same mapping built
  from different bases of the same kernel compared unequal - a trap, since
  nothing about a notation's meaning depends on which basis its kernel happens
  to be stored as. Struct equality was never actually what any caller wanted;
  the places that compared two notations for being "the same" now compare
  `mapping()` directly, which is what decides everything a notation does.
- **Accidental symbols are keyed to the prime, with one exception.** A
  notation keeping more than one accidental gives each its own fixed symbol
  from `PRIME_SYMBOLS`, so the same prime prints the same way regardless of
  what else is in the notation. A notation keeping exactly one is read as ups
  and downs does, for a generic step, so it gets the generic `^`/`v` - unless
  that one accidental is the quartertone `33/32`, which gets its own `t`/`d`
  instead. `accidental_symbol` in `notation.rs` is the whole of this.
- **Six symbols means 19-limit is the ceiling** for `from_ji`, and for any
  notation keeping more than one accidental generally. A single accidental
  never needs a symbol of its own - see above - so a notation with only one
  works at any prime.
- ~~`Temperament::et` silently defactors.~~ Fixed: `Temperament::from_mapping`
  now refuses a contorted mapping instead of silently saturating it, since
  saturating changes which temperament it describes. `et(24, 2.3.5)` is the
  case that showed it, contorted and saturating to `12et`; it now errors
  instead of returning `12et` for a question about `24et`. Checked by
  comparing `hnf(mapping)` to `saturation(mapping)`, which agree exactly when
  there is no contorsion. `from_commas` never needed the check: a kernel is
  always saturated.
- **An accidental defined as one step** - what ups-and-downs uses in general -
  is still not available, so an equal temperament with no `(2, 3, p)` accidental
  worth one step gets no notation rather than an arbitrary one. Over `2.3.5`
  that is 25, 51 and 54 up to 72; over `2.3.5.7`, 25, 54 and 57. A wider
  subgroup is the answer for these, and the `24et` entry shows it working.
- **`Weighting` is not exposed on `options`.** The fallback in step 2 hardcodes
  the default. Tenney and Wilson agree exactly - same commas, same spellings -
  over the temperament list and over every equal temperament to 99 in the 5, 7
  and 11-limit, so nothing turns on it yet. It is no longer likely to go away:
  two metric-free replacements were tried in the notation basis and both are
  worse, for the reason given there.
- **`shortest_stack` walks a box, like `Simplifier` does.** Minimising the
  number of marks is an `L1` problem again, so the same pattern applies: reduce
  the relations with `lll`, then walk `WIDTH` either way and settle ties on the
  counts. Nothing over the list or any equal temperament to 99 comes within the
  box of its edge - checked by instrumenting the walk - but it is a bounded
  search, not a proof, and a subgroup wide enough to need a larger `WIDTH` would
  fail quietly rather than loudly.
- **Whether `225/224` may be an accidental: resolved, no.** The `{2, 3, p}`
  support rule is a real constraint, not just an accident of how accidentals
  are derived. Everything hangs off a prime having exactly one accidental -
  `Search`'s necessary/optional/passed-over classification, the comma
  substitution being triangular, symbols keyed to a prime - and that structure
  is worth more than what a multi-prime accidental like `225/224` would buy.
  Nothing stops it existing in principle, but the search as built would never
  find one, and that is not worth changing.
