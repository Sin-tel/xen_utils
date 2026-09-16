# Notes

What the library does and why it does it that way. The code is the reference for
what is implemented; this is for the reasoning, and for what is still open.
Everything stated here was checked against the code at the time of writing.

`math.md` has the algebra. This does not repeat it.

## The model

A **temperament** is a surjection `T : Z^n -> Z^r` from the interval vectors of a
just intonation subgroup to pitches. It is the ground truth: it says which
intervals are the same note.

A **notation** is a legible interface to it. Notation coordinates are counts of
generators - the octave `2/1`, the fifth `3/2`, then one accidental per prime
worth keeping - and a written note is a vector of them. `F C G D A E B` are the
fifth coordinates `-1` to `5`; seven fifths beyond those is a sharp; each further
coordinate is an accidental mark.

Just intonation is a useful fiction about the temperament. The approximations
matter, but in 41et `14/11` and `81/64` **are** one note, and so are `E`, `vvvF`
and `vB#`.

## A notation is its generators

Everything a notation does follows from where its generators sit in the
temperament:

```
images      = temperament.map_all(generators)   // rank x r
enharmonics = kernel_left(images)               // rank - rank(temperament)
```

`images` is `Notation::pitch`, the map from a written note down to what it
sounds. Its kernel is the **enharmonic lattice**: the written notes the
temperament calls one pitch. `Notation::build` is those two lines plus two
checks - that every accidental beyond the first has a symbol, and that the
generators reach every pitch at all (`spans`).

**There is no map from just intonation to choose**, and that is the whole shape
of the design. Two notations keeping the same accidentals over the same
temperament are the same notation, so `generators()` is what to compare.

## The two questions

A pitch has two questions and they are independent. Neither constrains the other,
and the library answers them in different places.

**What is this note?** - ranked just intervals, `Simplifier::candidates`. Depends
on the temperament alone, so every notation of a temperament answers alike.

**How else can I write it?** - ranked spellings, `Notation::respell`. The
spellings of one pitch are a coset of the enharmonics, so this depends on the
generators and the temperament and on nothing else.

`Notation::spell` composes the second with the temperament: ask what pitch an
interval is, then ask which of that pitch's spellings reads best. In 41et
`spell(14/11)` and `spell(81/64)` therefore give the same answer, because they
are the same note.

`cargo run --example candidates` and `--example spellings` are the two lists.

## Accidentals

Derived, not given. For each prime beyond 3, shift it along the fifth chain by
0, 1, -1, 2, -2, ... fifths, octave-reduce each candidate to within a tritone of
a unison, and take the first that lands within **half an apotome**
(`MAX_ACCIDENTAL_CENTS`, 56.84 cents, `600 * (7 * log2(3) - 11)`).

That bound is load-bearing. Anything wider is closer to the neighbouring point
of the fifth chain, so it belongs there as a sharp or a flat instead. It is also
what makes a mark and a sharp commensurable further down.

```
5  -> 81/80     21.5c      13 -> 1053/1024  48.3c
7  -> 64/63     27.3c      17 -> 4131/4096  14.7c
11 -> 33/32     53.3c      19 -> 513/512     3.4c
```

The 13, 17 and 19 choices are the Helmholtz-Ellis ones, for free. `33/32` is the
widest in normal use, 3.6 cents inside the bound; `27/26` at 65.3 cents is what
the bound excludes for 13.

By construction an accidental has exponent `+-1` on its own prime and support
`{2, 3, p}`. Normalising them all to ascending gives the conventional directions
for free.

## Ranking spellings: a sharp is worth two marks

`apotomes` in `notation.rs` is the whole of it:

```
7 * marks + 2 * |fifths - 2|
```

Seven fifths are an apotome, so the chain and the accidental marks are
commensurable and nothing has to be tuned. An accidental is at most **half** an
apotome by construction, so two marks are an apotome, which is a sharp: seven
half-apotomes for a mark and two for a fifth.

The fifths are counted from `D`, the middle of the seven naturals, rather than
from `C`. Measuring from `C` makes the flat side cheaper than the sharp side by
one fifth throughout, which is enough to spell 5et's third `F` and 7et's `Eb`.
Ties fall back to the coordinates, so the order never depends on how the search
was walked. The octave does not appear: register costs nothing to write, and the
pitch fixes it once the rest is chosen.

Out of this, without being told anything about conventional practice:

```
12et   C  Db  D  Eb  E  F  F#  G  Ab  A  Bb  B
```

with the sharps as the second choice each time and the tritone an honest tie.

### What the other counts do

Three readings were tried and only the test suite separated them.

| measure | result |
| --- | --- |
| `7 * marks + \|f\|`, from `C` | 22et writes `D#`, never using the comma it keeps |
| `7 * marks + 2 * beyond the naturals` | differs on 12 of ~900 notations, all 11et and 13et, where it stops using an accidental the notation keeps |
| `7 * marks + 2 * \|f - 2\|` | in the code |

Two cheaper counts fail outright, and both fail the same way - by letting one
kind of symbol be preferred at any price:

- **Fewest marks, sharps only a tie-break.** An equal temperament whose fifth
  chain reaches every pitch always has, far out along it, a spelling with no
  marks and six sharps. Over 41et the first two entries survive a wider search on
  7 of 41 pitches.
- **Fewest symbols**, marks and sharps counted alike. Stable, but it cannot
  separate `C##` from `Ebb`, and it puts 41et's `B#` ahead of `^C`.

## The search: which accidentals to keep

That is all `Search` decides, and it is the only thing nothing downstream can do:
it picks the generators, hence the enharmonics, hence the notation. 226 lines.

Each accidental falls into one of four classes:

- **tempered out** - the temperament maps it to zero, so it would raise by
  nothing. Never kept.
- **passed over** - worth exactly what an accidental already kept is worth, up to
  direction, so it would be that accidental over again. `81/80` and `64/63` being
  one interval is a property of 41et, not a second symbol to read. Never kept.
- **necessary** - the smallest subset of what is left that reaches every pitch at
  all. Always kept. A lexicographic subset search, so the accidentals offered
  first are kept where there is a choice.
- **optional** - the rest. The run takes them on one at a time.

`Notation::options` is the run from keeping only the necessary ones to keeping
them all.

### The one place rank 1 is singled out

**An equal temperament's finest accidental must be worth one step.** That is the
rule ups and downs is built on: with no symbol for a single step, single steps
can only be reached by walking the fifth chain, which is not how anyone writes
such a temperament. So if no accidental is worth a single step, an equal
temperament is offered none at all.

It has to be singled out, because above rank 1 no accidental can reach every
pitch by itself and there is nothing for "worth one step" to generalise to.

The cost is that some equal temperaments get no notation rather than a bad one:
25, 51 and 54 over `2.3.5`; 25, 54 and 57 over `2.3.5.7`; 54 and 57 over
`2.3.5.7.11`. A wider subgroup is the answer - 24et is contorted over `2.3.5` but
over `2.3.5.11` its quartertone is `33/32`, worth exactly one step.

## Which notation to recommend

`from_temperament` returns **the smallest notation that keeps every nominal**:
one where each prime is written on the letter just intonation gives it. Each
further accidental is another symbol to read, so smaller is better - but not at
the price of moving a prime onto another letter.

41et is the case that makes it obvious:

```
   [2]                  5/4 Fb5   7/4 Cbb6   11/8 Abbb5
-> [3]  81/80           5/4 vE5   7/4 vBb5   11/8 ^^F5
   [4]  81/80 33/32     5/4 vE5   7/4 tA5    11/8 tF5
```

`[2]` has all three primes off their nominals. `[4]` only turns two marks into
one of another kind, so `[3]` is the one wanted.

**Sharps and flats do not count**, since seven fifths leave the letter alone.
Flattone writes `11/8` as `F#5` where just intonation has `^F5`, and that is
still an `F`, so there is no reason to take on an accidental for 11.

Where nothing qualifies the recommendation falls back to the last notation.

**This rule is the weakest thing in the library.** Nothing stops a notation
keeping its nominals by stacking a pile of marks, and nothing says the smallest
notation that manages it is the one anyone wants. It happens not to go wrong on
anything in the list, which is not the same as being right.

## Simplifying

Which interval a tempered pitch "is" is a closest vector problem on the
temperament's comma lattice. Three things about it took working out and are
settled.

**Which lattice.** The whole comma lattice, `ker T`. Not the enharmonics, which
stay inside what the notation's generators span and so could never answer `7/4`
in a notation over `{2, 3, 5}`.

**Which basis.** Just intonation coordinates. Wilson and Tenney weight *prime
exponents*, and there is no prime on the fifth axis, so weighting notation
coordinates would mean nothing.

**Which norm.** Wilson, `sopfr`, the sum of the prime factors of numerator and
denominator. It is an **L1** norm and the lattice algorithms take a quadratic
form, so the two disagree: under L2 a step of 41et is `45/44`, under L1 it is
`64/63`. So `lll` and `cvp_exact` seed the search and a bounded walk finishes it
under `sopfr`, one ball at a time out to `SEARCH_RADIUS`.

What comes out for 41et is the interval table anyone would write by hand:

```
1/1  64/63  25/24  21/20  16/15  12/11  10/9  9/8  8/7  7/6  32/27  6/5
11/9  5/4  14/11  9/7  21/16  4/3  27/20  11/8  7/5  10/7 ...
```

**The answer is the first of a list.** The walk passes over every interval in
every ball it steps through, so keeping them costs nothing and `candidates`
returns them sorted. That is a neighbourhood rather than the whole coset, so only
the right end of the list means anything.

One dead end, for the record: **leaving the octave out of the cost is wrong**.
Free powers of two let the search buy any number of them to avoid one high prime,
which inverts exactly the identifications wanted - 41et's `14/11` becomes
`81/64` and its `9/7` becomes `32/25`. `sopfr` is a statement about the ratio and
the 2s are part of the ratio.

## Performance

Built once, reused after, so constructing a `Notation` or a `Simplifier` is not a
concern. `cargo run --release --example bench` times what runs repeatedly, on
41et over `2.3.5.7.11`:

```
Notation::spell            2.4 us
Notation::to_just           39 ns
Simplifier::simplify        56 us
Simplifier::candidates(8)   54 us
```

`to_just` is a matrix multiply and not worth a second thought. **`spell` is no
longer free**: it used to be one too, at the same 39 nanoseconds, and it is now a
diophantine solve, a closest vector and a box walked under `apotomes`. Still
cheap enough to spell every visible note on a redraw - a hundred notes is a
quarter of a millisecond - but no longer in the same class, and not to be called
in a loop that does not need it.

It was 48 microseconds before the obvious fix, which is the same one
`Simplifier::search` needed once: ranking a candidate cloned its coordinates to
break the tie with.

## Everything here walks a box, and every box needs a reduced basis

`shortest_stack`, `respell`, `simplify` and the examples all do the same thing -
seed a point, then walk a bounded box around it under the norm actually wanted -
and the same two mistakes were made in each.

**Reduce the basis first, against the form the search measures with.**
Unreduced, 41et's enharmonics come back as "the octave is 41 ups" and "the fifth
is 24 ups", and a box around those finds `vvvvvvvD` at seven marks while never
reaching `vB#` at one mark and one sharp. Reducing against a *different* form
from the one the search uses is the subtler version of the same error: the basis
comes out short in the wrong sense. `Notation::weights` is built once and used
for both the reduction and the search, carrying the squares of what `apotomes`
costs - four for a fifth, forty nine for a mark.

**Seed near the answer.** `solve_diophantine` hands back any solution at all, and
it can be far out: a two-accidental notation once produced a spelling fifteen
marks from the best one, and orwell's replacement for `64/63` sat six relations
away from `225/224`.

`lll` reads `basis[0]` before it checks for an empty basis, so a notation that
spells bijectively has to skip the call.

## What the derivation used to do, and why it stopped

Worth keeping because it bounds what a search can usefully decide.

The notation used to be a linear map `N` from just intonation to notation
coordinates, built by stacking the generators on a fixed comma per dropped
accidental into a square matrix and inverting it. Choosing those commas was most
of the search: a preference order over stacks of the accidentals still kept, and
a metric fallback when the fifth chain had to be walked.

It was deleted because **nothing downstream could see which comma was picked.**
Two notations over 41et built from different comma rules have, checked by Hermite
normal form, *identical* enharmonic lattices - `[[1, 0, -41], [0, 1, -24]]` for
both - because `pitch` is fixed by where the generators go and the commas do not
enter it. Their kernels genuinely differ, so they differ in which just intervals
they spell alike; but every spelling either could give lies in the same coset of
the same lattice, and `respell` ranks that coset directly.

Three things were measured on the way and are worth not re-deriving:

- **The run never depended on the commas.** Comparing rank and accidentals kept
  across ~900 notations under two different comma rules: zero differing lines.
- **The preference order decided almost nothing.** Replacing it with any
  solution at all moved 17 of 912 lines, every one the same shape: the `33/32`
  comma of an 11-limit notation keeping both `81/80` and `64/63`.
- **Picking the spelling from the coset instead makes the comma irrelevant**, and
  is cheaper on the page - 28 marks against 33 over 41et. That is the design now.

The thing that had to be accepted to get here is that spelling stops being a
homomorphism. It costs nothing in practice because notation coordinates are the
state: transposing is a notation-coordinate operation and stays exact, and the
ranking only runs when someone asks. But spelling a chord note by note from just
intonation ranks each note independently, so anything wanting a consistent chord
should transpose in notation coordinates rather than respell each member.

## Loose ends and known limits

- **The ranking is blind to which accidental it uses.** 41et's rank 4 notation
  writes `7/4` as `tA5`, using the mark that belongs to 11, because one mark on
  `A` is cheaper than one on `Bb`. Flattone's larger notation keeps an accidental
  for 11 and writes `11/8` as `F#` anyway - which is right, since `F#` and `tF`
  are the same note and `F#` is simpler. Where the cheapest symbol carries the
  wrong harmonic hint this is taste, and there is no rule underneath it: `spell`
  is handed a pitch, not a prime, so it cannot prefer the accidental belonging to
  what is being written.
- **The recommendation is the weakest rule here.** See above.
- **Making `keeps_nominals` follow the accidental is vacuous**, though it looks
  like the fix. Asking for the prime on the letter the chain gives *plus*
  wherever this notation writes the accidental itself is a condition `spell`
  satisfies almost automatically, and 41et then recommends its bottom notation.
- **The simplifier is a local search, not a proof.** `SEARCH_RADIUS` is 2, and
  the note on it records that a radius of 2 walks as far as a radius of 3 for
  every interval within nineteen generators of the unison. A radius of 1 is
  provably not enough: `a_radius_of_two_is_enough_and_a_radius_of_one_is_not`
  pins the case, 9et reaching `1/39366` of norm 29 where a radius of 1 stops at
  `1/34560` of norm 30.
- **`respell` is a bounded search too.** `verify` checks that `spell` returns the
  cheapest spelling in its coset, walking wider than `respell` does and writing
  the cost out a second time so that it checks the answer rather than restating
  how it was found. 0 failures over the sweep, which is what says `RESPELL_WIDTH`
  and the seed are adequate.
- **Six symbols means 19-limit is the ceiling** for any notation keeping more
  than one accidental. A notation keeping exactly one reads it as ups and downs
  does and gets the generic `^`/`v`, so it works at any prime - unless that one
  accidental is `33/32`, which gets `t`/`d`.
- **An accidental defined as one step** - what ups and downs uses in general - is
  still not available, which is why some equal temperaments get no notation.
- **Whether `225/224` may be an accidental: no.** The `{2, 3, p}` support rule is
  a real constraint. Everything hangs off a prime having exactly one accidental -
  the necessary/optional/passed-over classification, symbols keyed to a prime -
  and that is worth more than a multi-prime accidental would buy.
- **An equal temperament sweep is a coverage test, not a judgement.**
  Temperaments can be arbitrarily bad and most of what a sweep turns up is
  nobody's notation. It catches panics and shows what a rule change moved; for
  whether an answer is sensible, read the named list, and the rank 2 and rank 3
  entries especially.

## Where the code lives

- `notation.rs` - the `Notation` type: `generators`, `pitch`, `enharmonics`,
  `spell`, `respell`, `note`, the ranking `apotomes`, and the derivation of a
  single accidental.
- `search.rs` - `Notation::options`, and only the choice of which accidentals to
  keep. Tested through `options`.
- `simplify.rs` - `Simplifier`: the comma lattice reduced once, then a seeded
  walk per interval.
- `temperament.rs`, `primes.rs`, `util.rs` - the temperament, the subgroup, and
  integer vector helpers with no music in them.

Examples: `verify` sweeps the invariants and prints a failure count; `dump`
prints a fingerprint of every notation to be diffed across a change; `notations`
lays out each run; `candidates` and `spellings` are the two questions;
`simplify` and `miracle` walk one temperament; `bench` times the hot path.
