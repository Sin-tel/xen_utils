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

### A dead end, recorded so it is not repeated

An earlier attempt reduced each dropped accidental independently to the
cheapest enharmonically equal spelling, with a weighted cost (fifths cost 1,
accidentals cost `w`). The weight traded two things off against each other:

| accidental weight | 2 | 3 | 4 | 5 | 6 | 7 | 10 |
|---|---|---|---|---|---|---|---|
| `ker` nesting failures | 2 | 2 | 2 | 2 | 0 | 0 | 0 |
| 41et rank 3, `11/8` | `^^F` | `^^F` | `^^F` | `^^F` | `vGb` | `vGb` | `vGb` |

Low weights gave the conventional `^^F` but broke the property that a smaller
notation's kernel contains a larger one's, so a rank 3 spelling could not be
expected to simplify onto the rank 2 one. High weights nested but spelled
`11/8` as `vGb`. **The cost function was the wrong tool** - the equal
temperament heuristic below gets `^^F` with no weights at all. Prefer rules
that are forced over rules that are tuned.

## Equal temperaments (implemented)

Two questions settle everything, with no costs or weights anywhere.

1. **Does the fifth chain reach every note by itself?** Yes exactly when
   `gcd(divisions, fifth_steps) == 1` - the circle of fifths comes back round
   only after all of them. If so, the fifth chain alone is already a notation.
2. **Is any accidental worth a single step?** Scan them in order of prime and
   take the first. An equal temperament notated with one accidental wants that
   accidental to be the step, so this is the one to reach for before any other.

From those:

- chain reaches everything -> the step accidental is **optional**
- chain does not -> the step accidental is **necessary**
- no step accidental and the chain does not reach -> **no notation**, bail

Then classify the rest by step count: **0** is tempered out and would raise by
nothing; **±1** is the step over again (41et has `81/80` and `64/63` both worth
one, and either serves, giving literally the same mapping); **more than 1** is
a further accidental worth having, offered in order of prime.

Replacement needs no choice because the step is worth exactly one step:

- with the step kept, a dropped accidental worth `k` steps becomes `k` of it
- with nothing kept, it becomes the interval on the fifth chain nearest a
  unison that the temperament cannot tell it from - one axis, so "nearest" is
  unambiguous

41et in the 11-limit therefore offers exactly three notations:

```
-               [2] 5/4 Fb5  7/4 Cbb6  11/8 Abbb5     the fifth chain alone
81/80           [3] 5/4 vE5  7/4 vBb5  11/8 ^^F5      ups and downs
81/80 33/32     [4] 5/4 vE5  7/4 vBb5  11/8 >F5       plus one for 11
```

`Notation::options` returns the run, smallest first; `from_temperament` picks
an end of it.

## Rank 2 and up: where it stands and how to move forward

Currently rank 2 and up only does the `r == rank(temperament)` subset search
(forced, no heuristic) and otherwise falls back to the notation keeping every
accidental the temperament does not temper out. So the two ends exist but the
middle does not, which only shows up for the three list entries with two
optional accidentals: **huygens**, **11-limit miracle** and **11-limit marvel**.
Those are the test cases.

Ideas, roughly in order of how much they buy:

1. **Transplant the equal temperament structure.** The same three classes
   apply: *necessary* is the kept set from the span search, *optional* is the
   dropped accidentals the temperament does not temper out (in order of prime),
   *useless* is the ones it tempers out. `options` then runs from enabling none
   to enabling all. The two ends already agree with what is built today, so
   this is mostly bookkeeping plus a replacement rule for the middle.

2. **Generalise "k copies of the step" to "a stack of the kept accidentals".**
   The equal temperament rule works because the step is worth one, so
   `T(a_j) = k * T(step)` has a unique solution. The general version is
   `T(a_j) = sum m_i T(a_i)` over the kept accidentals, unique whenever their
   images are independent. When it solves, the prime keeps the nominal and
   sharps that just intonation gives it and only the accidental count changes.
   No weights. This is exactly the equal temperament rule when
   `rank(temperament) = 1`, since a step accidental's image generates all of
   `Z` and so every image lies in it.

   Where it does not solve - 11-limit miracle keeping `81/80` cannot express
   `64/63`, and vice versa - **fall back to the next larger notation** rather
   than inventing a spelling. Refusing to guess is better than guessing.

3. **Fix an order and construct top-down.** Give each accidental one comma
   whose replacement may use the octave, the fifth, and accidentals *later* in
   a fixed order (highest prime first). Then:
   - kernel nesting is automatic, since disabling more accidentals only adds
     commas, so a smaller notation is always a coarsening of a larger one;
   - unimodularity is automatic for *any* subset, because substituting
     `comma_j = a_j - c_j` is triangular in the order and terminates at the
     octave and fifth, recovering the just intonation generators;
   - each replacement is computed with the largest set of accidentals
     available, which is what makes `33/32` become `2 x 81/80` rather than
     something far out on the fifth chain.

4. **Leave prettiness to `simplify`.** Any notation above the temperament's
   rank has a nonzero enharmonic lattice, and the mapping only needs *some*
   representative. Do not spend effort on the raw output beyond the forced and
   weight-free rules; reducing a spelling to the nicest of its enharmonic
   equivalents is a separate feature that will need a cost function of its own
   (and that is the right place for one).

## Loose ends and known limits

- **`simplify` is not started.** The pieces it needs are in place: the
  enharmonic lattice is `kernel_left` of the generator images, and reduction is
  a small closest-vector problem in a metric that ignores the octave
  coordinate. `diophantine` has `lll`, `nearest_plane` and `cvp_exact`; a
  bounded box search after reduction also works and takes an arbitrary cost
  function, which matters because the natural cost is not quadratic (an
  asymmetric nominal window, accidentals weighted against sharps).
- **Derived `PartialEq` on `Notation` compares the stored comma basis**, so two
  notations with the same mapping built from different bases of the same kernel
  compare unequal. `options` works around this by comparing mappings. Either
  implement `PartialEq` by hand or keep remembering.
- **Accidental symbols are positional**, so the same comma prints as `^` in one
  notation of a temperament and `>` in another. Fine for debugging output;
  wrong if symbols should be stable across a temperament's run, in which case
  key them to the prime.
- **Four symbols means 13-limit is the ceiling** for `from_ji`. Nothing else
  limits the subgroup.
- **`Temperament::et` silently defactors.** `et(24, 2.3.5)` saturates to
  `[12, 19, 28]`, i.e. to 12et, because the 24et patent mapping of the 5-limit
  is contorted. 24et as such cannot be asked for this way.
- **Equal temperaments with no notation**: over `2.3.5` these are 25, 51 and 54
  up to 72; over `2.3.5.7`, 25, 54 and 57. Their fifth chains close early and
  no `(2, 3, p)` accidental is worth one step - 25et closes after five notes
  and its syntonic comma is worth two steps. Handling them needs an accidental
  that is defined as one step rather than as a comma, which is what
  ups-and-downs does in general and is out of scope for now.
- **Whether `225/224` may be an accidental** was left open. Nothing enforces
  the `{2, 3, p}` support rule, because accidentals are derived rather than
  given, and the derivation cannot produce anything else.
