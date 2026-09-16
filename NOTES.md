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

## A notation is its generators

The data is a subgroup, a temperament and a list of accidentals. Everything
else is derived from where the generators - the octave, the fifth and those
accidentals - sit in the temperament:

```
images[i]   = temperament.map(generator[i])     // rank x r
enharmonics = kernel_left(images)               // rank r - rank(temperament)
```

`images` is the map `pitch` from a written note down to what it sounds, and its
kernel is the enharmonic lattice: the written notes the temperament calls one
pitch. `Notation::build` is those four lines plus two checks - that every
accidental beyond the first has a symbol, and that the generators reach every
pitch at all, which is `spans`.

**There is no map from just intonation to choose.** Spelling a just interval is
two steps that were previously one: ask the temperament what pitch it is, then
ask which of that pitch's spellings reads best. That is `spell`, and the
alternatives it passed over are `respell`. Both depend on the generators and the
temperament and on nothing else, so **two notations keeping the same accidentals
are the same notation** - which is why `generators()` is what to compare and why
the search has nothing to say about just intonation.

This replaces a matrix inversion, a kernel spanned by one fixed comma per
dropped accidental, and about a hundred and ten lines of search picking those
commas. See "What a notation decides, and what it cannot" for why none of it was
deciding anything a caller could see.

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

## What a notation decides, and what it cannot

Worth settling, because it bounds what the search is for.

**Compare lattices by Hermite normal form, never by a printed basis.** A basis is
not canonical and two different ones span the same lattice. Everything below is
an HNF comparison.

### Two notations with the same generators share their enharmonics

Take 41et over `2.3.5.7.11` and its rank 3 notation, derived once as the code
does it and once with step 1 disabled so the simplest comma runs everywhere.
Both keep the octave, the fifth and `81/80`. Then the map `t` from notation
coordinates down to pitches is determined by where those three generators go,
which is the same for both - so `E = ker t` is the same lattice, and it is:

```
[1, 0, -41]      the octave is 41 ups
[0, 1, -24]      the fifth is 24 ups
```

identical for both, in HNF. The same holds for 31et. So as **symbol systems** -
what marks exist, what each is worth, which written notes name one pitch - the
two are the same, and any two spellings of a pitch differ by an element of `E`
whichever notation produced them.

### But their kernels differ, so they are not the same map

```
derived:          [1, 15, -5, -1, -3]      [0, 52, -17, -3, -10]
simplest comma:   [1,  2, -6,  1,  2]      [0,  6, -14,  2,   5]
```

Different lattices in HNF. Since `ker N + E = ker T` and `E` is shared, the
difference is entirely in `ker N` - **which just intervals the notation spells
alike**. A notation is not only a set of symbols; it is a map from intervals to
them, and that map is what the commas choose.

### The pitch question and the interval question

These are different questions and only one of them is the notation's.

- **Asked for a pitch**, the notation decides nothing. 12et's second step may be
  `D`, `Ebb` or `C##` and nothing in the notation prefers one. That is a coset of
  `E` and it needs a rule of its own.
- **Asked for an interval**, the notation decides everything, and the answer
  carries the harmonic reading the temperament discarded. In 12et, still with no
  freedom at all in its kernel:

```
16/15  -> Db        25/24  -> C#        135/128 -> C#
 9/8   -> D         10/9   -> D
```

  One step, two spellings, and the notation tells them apart because it was
  asked about intervals. `cargo run --example spell` is that table.

The pitch question is already answered upstream, by composing the two: 12et's
second step simplifies to `9/8` and so is written `D`; its first simplifies to
`16/15` and so is written `Db`, not `C#`. That is the right place for it.

### Choosing the spelling from the coset instead does not work

If the spelling were picked from the coset of `E` by counting symbols, the comma
choice could not affect it at all - both notations offer the same coset. It is
cheaper on the page: over 41et it costs 28 marks against the derived notation's
33. And it is wrong, because it answers the pitch question while the simplifier
is answering the interval question, and the two then disagree:

```
step 14   simplest 14/11   derived vdF   coset E
step 27   simplest 11/7    derived ^tG   coset Ab
```

A reading of `14/11` displayed as a plain `E` is a pythagorean third on the page
and an undecimal one in the analysis. The spelling has to be a function of the
interval, and that function is `N`. `cargo run --example table` prints both
columns.

### So what the search is for, and what it now is

**Implemented.** `search.rs` went from 470 lines to 220 and no longer mentions a
comma. What is left is: derive one accidental per prime, drop the ones the
temperament tempers out, drop the ones worth what a kept one is worth, and check
the rest one at a time against `spans`. `Plan` is two lists rather than four.


- **Which accidentals to keep** - tempered out, passed over, necessary,
  optional. This picks the generators, hence `t`, hence `E`: it chooses the
  symbol system itself. Nothing downstream can do this and it is the bulk of
  `Search`.
- **The run and the recommendation**, which are built on that.
- ~~The commas~~, which used to pick one map out of the many compatible with
  that system. Gone: nothing downstream could see which one was picked, since
  every spelling it could have given is in the same coset of `E` and `respell`
  ranks that coset directly.

### The question that is actually wanted

Settled: **given a spelling, what else does the temperament call the same
note?** The kernel of the notation does not enter it. The spellings of one pitch
are a coset of `E`, and `E` is `kernel_left` of the generator images - so it
depends on the generators and the temperament and on nothing else, which is why
two notations built on the same generators answer identically and why the
simplifier does not care. `Notation::enharmonics` already computes it that way.

The temperament is the universe and the notation is a legible interface to it.
Just intonation is a useful fiction about that universe: the approximations
matter, but `14/11` and `81/64` are one note in 41et, and so are `E`, `vvvF`,
`^^D#` and `vB#`. Two separate queries, then, and both are wanted:

- **what is this note?** - ranked just intervals, `Simplifier::candidates`,
  which depends on the temperament alone;
- **how else can I write it?** - ranked spellings, the coset of `E`, which
  `cargo run --example respell` prototypes.

Which two of four commas a rank 3 notation of 41et should lose is not a question
to put to anyone, and it does not have to be asked: it only decides which
spelling comes back first for a just interval, and that is what the second query
is for.

### Ranking spellings: a sharp is worth two accidental marks

Seven fifths are an apotome, so the chain of fifths and the accidental marks are
commensurable and nothing has to be tuned. And an accidental is at most **half**
an apotome by construction - `MAX_ACCIDENTAL_CENTS`, the rule the whole
derivation of an accidental hangs on - so two marks are an apotome, which is a
sharp. In half-apotomes: **seven for a mark, two for a fifth**.

The fifths are counted as how far the note sits *beyond the seven naturals*
rather than as a distance from `C`. `F C G D A E B` are the fifth coordinates
`-1` to `5`, and all seven cost nothing; a fifth past either end costs two.
Measuring from `C` instead makes the flat side cheaper throughout, which is
enough to spell 5et's third `F` and 7et's `Eb`. Ties are settled by distance
from `D`, the middle of the naturals, and then by the coordinates, so the order
never depends on how the search was walked.

That is `apotomes` in `notation.rs`, and it is the whole of the ranking. It
gives 12et `C Db D Eb E F F# G Ab A Bb B` with the sharps as the alternates,
41et's one step `^C` ahead of `B#`, and 22et's third `vE` ahead of `D#` in the
notation that keeps a comma for it.

Three readings of "a fifth is a seventh of an accidental" were tried against the
test suite, which is the only thing that separated them:

| measure | what it costs a fifth | result |
| --- | --- | --- |
| `7 * marks + \|f\|` | one, from `C` | 22et writes `D#`, never using the comma it keeps |
| `7 * marks + 2 * \|f - 2\|` | two, from `D` | flattone writes `F#`, 7et `Eb` |
| `7 * marks + 2 * beyond` | two, past the naturals | what is in the code |

### The gap this leaves

`Notation::respell` now answers it: given a written note, the other ways to
write that pitch, best first. `candidates` cycles *readings* - distinct just
intervals - and `respell` cycles spellings, which are the two questions a caller
has and the two the library now separates.

## Which notation to recommend

`from_temperament` no longer takes an end of the run. It returns **the smallest
notation that keeps every nominal**: one where each prime is spelled on the
letter just intonation gives it. Each further accidental is another symbol to
read, so smaller is better - but not at the price of moving a prime onto another
letter.

41et is the case that makes it obvious once the run is laid out:

```
   [2]  5/4 Fb5   7/4 Cbb6   11/8 Abbb5
-> [3]  5/4 vE5   7/4 vBb5   11/8 ^^F5
   [4]  5/4 vE5   7/4 tA5    11/8 tF5
```

`[2]` has all three primes off their nominals; `[4]` only turns the two marks of
`[3]` into one of another kind, so `[3]` is the one wanted.

Since `keeps_nominals` now asks where `spell` puts a prime rather than where a
comma did, it can answer differently: huygens' rank 3 notation already spells
every prime on its own nominal, so the recommendation moved down one.

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

## Where the code lives

- `notation.rs` - the `Notation` type: its generators, `pitch`, `enharmonics`,
  `spell`, `respell`, `note`, the ranking `apotomes`, and the derivation of a
  single accidental.
- `search.rs` - everything behind `Notation::options`, and only the choice of
  which accidentals to keep. `Search` holds the temperament, the accidentals and
  their images; `Plan` sorts them into necessary and optional. Tested through
  `options`.
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

- **`respell` walks a box, like everything else here.** `RESPELL_WIDTH` either
  way around a `cvp_exact` seed on the `lll` reduced enharmonics. Both the
  reduction and the seed are load-bearing and both were found the hard way:
  unreduced, 41et's enharmonics are "the octave is 41 ups" and "the fifth is 24
  ups" and a box around them misses everything legible; unseeded, a notation
  with two accidentals hands back a spelling fifteen marks out. The quadratic
  weights `cvp_exact` takes are a stand-in for `apotomes`, which is a count and
  not a form.
- **The ranking is blind to which accidental it uses**, and that shows. 41et's
  rank 4 notation writes `7/4` as `tA5`, using the mark that belongs to 11,
  because one mark on `A` is cheaper than one mark on `Bb` - `A` is inside the
  naturals and `Bb` is a fifth outside. Flattone's larger notation keeps an
  accidental for 11 and then writes `11/8` as `F#` anyway. Both follow from
  spelling being a function of the pitch alone, which is what was wanted, and
  both are cases where the cheapest symbol carries the wrong harmonic hint. A
  tie-break preferring the accidental of the prime being written is not
  available, since `spell` is handed a pitch and not a prime; a cost that reads
  the accidentals in order of their primes would be.
- **Huygens now wants the second of its three notations**, not the third, its
  rank 3 already spelling every prime on its own nominal. That is
  `keeps_nominals` asking `spell` rather than asking a comma.
- ~~Nothing offers the other spellings of a pitch.~~ Fixed: `Notation::respell`.

- **The simplifier is a local search, not a proof.** `Simplifier` reduces the
  temperament's comma lattice once, takes the `cvp_exact` answer as a seed, and
  then walks it under the norm actually wanted. It agrees with a wider step for
  every interval within nineteen generators of the unison, over 59109 swept; the
  stalls past that are twenty octaves out, where the norm is minimized by
  trading a pile of one prime for a pile of another. See `SEARCH_RADIUS`.
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
- **Whether `225/224` may be an accidental: resolved, no.** The `{2, 3, p}`
  support rule is a real constraint, not just an accident of how accidentals
  are derived. Everything hangs off a prime having exactly one accidental -
  `Search`'s necessary/optional/passed-over classification, the comma
  substitution being triangular, symbols keyed to a prime - and that structure
  is worth more than what a multi-prime accidental like `225/224` would buy.
  Nothing stops it existing in principle, but the search as built would never
  find one, and that is not worth changing.
