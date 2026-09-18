# Notes

What the library does and why it does it that way. The code is the reference for
what is implemented; this is for the reasoning, and for what is still open.
Everything stated here was checked against the code at the time of writing.

## The model

A **temperament** is a surjection `T : Z^n -> Z^r` from the interval vectors of a
just intonation subgroup to tempered intervals. It is the ground truth: it says
which intervals are the same note.

A **notation** is a legible interface to it. Notation coordinates are counts of
generators - the octave `2/1`, the fifth `3/2`, then one accidental per prime
worth keeping - and a written note is a vector of them. `F C G D A E B` are the
fifth coordinates `-1` to `5`; seven fifths beyond those is a sharp; each further
coordinate is an accidental mark.

Just intonation is a useful fiction about the temperament. The approximations
matter, but in 41et `14/11` and `81/64` **are** one note, and so are `E`, `vvvF`
and `^^D#`.

## Words

Three spaces, one name each, and each vector is only ever called by its name:

- an **interval** is a just interval: prime exponents over the subgroup, of
  length `Subgroup::dim`;
- a **spelling** is notation coordinates: the octave, the fifth, then each
  accidental, of length `Notation::len`. `note` prints one;
- a **tempered** interval is generator counts in the temperament, of length
  `Temperament::rank`.

**Pitch** is kept for a real size in cents, which needs a tuning the library
does not have yet.

Mapping **down** is `temper`, from an interval or a spelling, and needs no
choice. Going **up** chooses an element of a coset, always starts from a tempered
interval, and is named for what it chooses: `spell` and `spellings`, `simplify`
and `simplifications`. Each has an `_interval` form that tempers a just interval
first. `Temperament::preimage` is an arbitrary element of the coset, not a
choice, and `Notation::to_interval` is a spelling's literal reading, taking each
accidental as the ratio it is named for.

## A notation is its generators

Everything a notation does follows from where its generators sit in the
temperament:

```
images      = temperament.temper_all(generators)   // len x r
enharmonics = kernel_left(images)                  // len - rank(temperament)
```

`images` is `Notation::temper`, the map from a spelling down to what it sounds.
Its kernel is the **enharmonic lattice**: the spellings the temperament makes
the same. `from_accidentals` is those two lines plus two checks - that no two
symbols collide, and that the generators reach every tempered interval at all
(`spans`).

**There is no map from just intonation to choose**, and that is the whole shape
of the design. Two notations keeping the same accidentals over the same
temperament are the same notation, so `generators()` is what to compare.

## The two questions

A tempered interval has two questions and they are independent. Neither
constrains the other, and the library answers them in different places.

**What is this note?** - ranked just intervals, `Simplifier::simplifications`.
Depends on the temperament alone, so every notation of a temperament answers
alike.

**How can I write it?** - ranked spellings, `Notation::spellings`. The spellings
of one tempered interval are a coset of the enharmonics, so this depends on the
generators and the temperament and on nothing else.

`spell_interval` composes the second with the temperament: temper the interval,
then ask which of its spellings reads best. In 41et `14/11` and `81/64` therefore
give the same answer, because they are the same note.

`cargo run --example simplifications` and `--example spellings` are the two
lists.

## Accidentals

Derived, not given. For each prime beyond 3, shift it along the fifth chain by
0, 1, -1, 2, -2, ... fifths, octave-reduce each candidate, and take the first that
lands within half a sharp.
(`MAX_ACCIDENTAL_CENTS`, 56.84 cents, `1200*log2(sqrt(2187/2048))`).

```
5  -> 81/80      21.5c
7  -> 64/63      27.3c
11 -> 33/32      53.3c
13 -> 1053/1024  48.3c
```

The 13, 17 and 19 choices are the Helmholtz-Ellis ones, for free. `33/32` is the
widest in normal use, 3.6 cents inside the bound; `27/26` at 65.3 cents is what
the bound excludes for 13.

By construction an accidental has exponent `+-1` on its own prime and support
`{2, 3, p}`. Normalising them all to ascending gives the conventional directions
for free.

## Ranking spellings: a sharp is worth two marks

`spelling_cost` in `notation.rs` is the whole of it:

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
tempered interval fixes it once the rest is chosen.

## Nominals are degrees

Give a written note a **degree**, the number of letters it sits up:
`7 * octaves + 4 * fifths`, and nothing for an accidental. A sharp is seven
fifths less four octaves, degree zero, so it keeps the letter too.

Just intonation gives each prime a degree of its own, and together those are a
linear map from interval vectors to degrees, `D = (7, 11, 16, 20, 24, 26, ..)`.
Its kernel holds `2187/2048` and every derived accidental, which is exactly why
those keep the letter.
`just_nominals` has it, with what the just spelling costs.

A prime is **on its nominal** when it is written at degree exactly `D(p)`.
Exactly, not modulo seven: in 41et `E##########` has the letter of `5/4` but
sits an octave below it.

`nominal_spellings` is the cheapest spelling of each prime at its just degree,
or `None` where no spelling of the tempered prime has that degree, and
`nominal_costs` what those cost. Spellings at a fixed degree differ by the
degree-zero enharmonics, so it is its own closest vector search on that slice -
it does not read what `spell` happened to pick.

A tempered interval can be written at `d + g Z`, `g` the gcd of the enharmonics' degrees.
Where `g > 1` (41et's fifth chain, `g = 4`) or no enharmonic moves the degree
(schismatic's), some letters are out of reach and the question has teeth. Where
`g = 1` every letter is reachable and only the cost tells notations apart.

**`keeps_nominals`**: every prime can be written on its nominal at a cost no
worse than either just intonation spends on it or this notation spends on its
cheapest spelling of it. The first half lets 41et's largest notation keep `vBb`
though `tA` is cheaper; the second lets 41et's middle one write `11/8` as `^^F`,
two marks where just intonation has one, because it has nothing better.

Only an accidental's image affects any of this. The ratio names it and is what
`to_interval` reads back; in 41et `49/48` is as good a step as `81/80`, though
just intonation writes `49/48` on `D`, so `to_interval(^C)` is letter-inconsistent.

## The search: which accidentals to keep

That is all `NotationOptions` decides. Before searching it drops an accidental:

- **tempered out** - the temperament maps it to zero, so it would raise by
  nothing.
- **passed over** - worth exactly what an earlier accidental is worth, up to
  direction, so it would be that accidental over again. `81/80` and `64/63` being
  one interval is a property of 41et, not a second symbol to read.

Every subset of the rest is a candidate, if it reaches every tempered interval. Of
each size the best is the one with

1. the fewest primes failing `keeps_nominals`,
2. then the lowest total of `nominal_costs`, a prime that cannot reach its
   nominal at all counting as worse than any cost,
3. then those costs compared one prime at a time, lowest prime first,
4. then the subset first in the order the accidentals were given.

`Notation::options` is the best of each size, smallest first, and
`Notation::with_count` the best of one size. The best of one size need not keep
what the best of the size below keeps; on everything `verify` sweeps it does,
and `verify` prints it where not.

Miracle is where the ranking earns its keep. Its octave and fifth leave a
quotient of `Z/6`, six generators to the fifth, and `81/80` is one generator,
`64/63` two and `33/32` three: `64/63` alone reaches only half the tempered intervals. Every
single accidental leaves 7 and 11 off their nominals, and the pairs tie three
ways on a total of 46, each writing one of 5, 7 and 11 with two marks. The
lower primes win, so it is `81/80 64/63` and `11/8 = ^>F`.

The candidates are the derived accidentals, or a list supplied to
`options_with`. Searching wider was tried in a scratch example and turned up
nothing better on miracle or marvel; supplying the list is the way to ask.

### The one place rank 1 is singled out

**An equal temperament's finest accidental must be worth one step.** That is the
rule ups and downs is built on: with no symbol for a single step, single steps
can only be reached by walking the fifth chain, which is not how anyone writes
such a temperament. So if no accidental is worth a single step, an equal
temperament is offered none at all.

It has to be singled out, because above rank 1 no accidental can reach every
tempered interval by itself and there is nothing for "worth one step" to
generalise to.

The cost is that some equal temperaments get no notation rather than a bad one:
25, 51 and 54 over `2.3.5`; 25, 54 and 57 over `2.3.5.7`; 54 and 57 over
`2.3.5.7.11`. A wider subgroup is the answer - 24et over `2.3.5` is not
primitive, saturating to 12et, but over `2.3.5.11` its quartertone is `33/32`,
worth exactly one step.

## Which notation to recommend

`from_temperament` returns **the smallest notation that keeps every nominal**.
Each further accidental is another symbol to read, so smaller is better - but not
at the price of moving a prime onto another letter.

41et is the case that makes it obvious:

```
   [2]                  5/4 Fb5   7/4 Cbb6   11/8 Abbb5
-> [3]  81/80           5/4 vE5   7/4 vBb5   11/8 ^^F5
   [4]  81/80 33/32     5/4 vE5   7/4 tA5    11/8 tF5
```

`[2]` cannot reach any of the three nominals. `[4]` only turns two marks into
one of another kind, so `[3]` is the one wanted.

Where nothing qualifies the recommendation falls back to the first notation. On
the data files that no longer happens for any named temperament: semaphore,
diaschismic with 17 and pele each had a cheaper spelling off the nominal than on
it, and the old rule, reading only `spell`, took that as a miss. Swept against
the old rule, the new one agrees everywhere the old one found anything.

## Simplifying

Which just interval a tempered interval "is" is a closest vector problem on the
temperament's comma lattice.

**Which norm.** Wilson, `sopfr`, the sum of the prime factors of numerator and
denominator. It is an **L1** norm and the lattice algorithms take a quadratic
form, so the two disagree: under L2 a step of 41et is `45/44`, under L1 it is
`64/63`.

**The search is exact.** `sopfr` is L1 with the primes as weights, and the L2
norm under `diag(p²)` never exceeds it: `sqrt(sum (p e)²) <= sum p |e|`. So
once some interval of `sopfr` `S` is found, everything at least as simple lies
in the L2 ball of radius `S`, and a Schnorr–Euchner enumeration of the comma
lattice that prunes on the current `count`-th best finds the true top `count`.
That is `diophantine::cvp_l1_top_k`; ties go to the quadratic norm, then to the
comma. The first leaf is the L2 closest point, so the bound is finite at once,
and the gap it has to close is at most `sqrt(n)`.

It replaced a walk of balls of `5^rank` points around the L2 seed, repeated
until nothing improved. That was a local search - 9et once needed two steps at
once to reach `1/39366` from `1/34560` - and exponential in the lattice's rank:
11 milliseconds a call at the 19-limit, where the enumeration takes about 10
microseconds. On 20,674 tempered intervals up to the 19-limit the two agreed on
every best answer, and the walk's top 8 was worse in 1,046 lists, since it only
ranked what it happened to pass.

## Performance

`cargo bench` (criterion, `benches/notation.rs`) times what runs repeatedly,
averaged over every step of the temperament. Per call:

```
                                41et 2.3.5.7.11   72et 2.3.5.7.11   41et 2.3.5.7.11.13
Notation::spell                     3.3 us            3.2 us             2.9 us
Notation::spellings(4)              3.0 us            3.0 us             3.0 us
Notation::spell_interval            3.0 us            3.0 us             3.0 us
Notation::to_interval                38 ns             38 ns              39 ns
Notation::temper                     34 ns             34 ns              34 ns
Simplifier::simplify                4.0 us            4.1 us             5.5 us
Simplifier::simplifications(8)      8.6 us            8.3 us              16 us
Notation::from_temperament          142 us            504 us             178 us
```

`to_interval` and `temper` are matrix multiplies and not worth a second thought.
**`spell` is not free**: it is a diophantine solve, a closest vector and a box
walked under `spelling_cost`. Still cheap enough to spell every visible note on a
redraw - a hundred notes is a third of a millisecond - but not to be called in a
loop that does not need it. `spellings` costs the same, since the solve and the
closest vector are most of it.

Simplifying is an enumeration and grows with the lattice's rank and with
`count`, slowly: about 10 microseconds for the best at the 19-limit.

Building a notation is once per temperament, but `from_temperament` searches
every subset of the accidentals and scores each, so it grows with them: 72et
keeps three of three and is the slowest here at half a millisecond.

`spell` was 48 microseconds before two fixes. The first was the one the
simplifier's old walk had needed too: ranking a candidate cloned its coordinates
to break the tie with. The second was narrowing `RESPELL_WIDTH` from 5 to 2, which
is where `spellings` asked for more than one answer spends almost everything -
that call went from 8.4 microseconds to 1.3.

A special case for one answer, which is what `spell` asks for, was worth 3.4x at
width 5 and 13% at width 2, so it is gone: most callers want several anyway.

## Everything here walks a box, and every box needs a reduced basis

`spellings` and the examples do the same thing - seed a point, then walk a
bounded box around it under the norm actually wanted - as `simplify` did until it
became an exact enumeration, and the same two mistakes were made in each.

**Reduce the basis first, against the form the search measures with.**
Unreduced, 41et's enharmonics come back as "the octave is 41 ups" and "the fifth
is 24 ups", and a box around those finds `vvvvvvvD` at seven marks while never
reaching `vB#` at one mark and one sharp. Reducing against a *different* form
from the one the search uses is the subtler version of the same error: the basis
comes out short in the wrong sense. `Notation::weights` is built once and used
for both the reduction and the search, carrying the squares of what `spelling_cost`
costs - four for a fifth, forty nine for a mark, and **nothing for an octave**.
It once charged four for an octave too, which pulled the seed towards `C5` and
made the best spelling depend on register: 2169 of 309094 calls in the width
sweep got even the single best answer wrong at width 1. The form is
then only semidefinite, but no enharmonic is a stack of octaves, so it is
definite on the enharmonic lattice and LLL and the closest vector are unaffected.

**Seed near the answer.** `solve_diophantine` hands back any solution at all, and
it can be far out: a two-accidental notation once produced a spelling fifteen
marks from the best one, and orwell's replacement for `64/63` sat six relations
away from `225/224`.

## Loose ends and known limits

- **The ranking is blind to which accidental it uses.** 41et's rank 4 notation
  writes `7/4` as `tA5`, using the mark that belongs to 11, because one mark on
  `A` is cheaper than one on `Bb`. Flattone's larger notation keeps an accidental
  for 11 and writes `11/8` as `F#` anyway - which is right, since `F#` and `tF`
  are the same note and `F#` is simpler. Where the cheapest symbol carries the
  wrong harmonic hint this is taste, and there is no rule underneath it: `spell`
  is handed a tempered interval, not a prime, so it cannot prefer the accidental
  belonging to
  what is being written. `nominal_costs` is where intent shows: it asks for the
  prime on its letter, whatever `spell` would choose.
- **The recommendation is still a rule, not a proof.** It agrees with every
  earlier answer and has no weight in it, but it reads only the primes. An
  interval such as `7/5` could be off its nominal in a notation that keeps every
  prime on its own.
- **The simplifier ranks only what it is asked for.** `simplifications` is an
  exact top `count`, so the search closes only once it holds `count` intervals:
  asking for all of them, say `usize::MAX`, never returns.
- **`spellings` is a bounded search.** The seed is an exact closest vector,
  so the walk only corrects for the quadratic form standing in for
  `spelling_cost`. `cargo run --release --example spellings_width` measures how
  far that correction reaches, against a much wider walk, over every notation
  of both data files and the equal temperaments to 72 (154547 calls):

  ```
  width   top 1 wrong   top 3 wrong   top 5 wrong
    0         13628        all           all
    1            35       1952         67643
    2             0          0           961
  ```

  Width 1 misses real answers, not ties - pele writes a tempered interval `vdF##`
  where `E##` is cheaper - so `RESPELL_WIDTH` is 2, pinned by
  `spellings_walk_wide_enough`. That is exact for the first three answers on
  everything swept, and approximate beyond: no fixed width can be exact for
  every `count`. An exact version would bound the walk by cost instead - the L2
  form never exceeds `spelling_cost`, so every spelling within cost `C` lies in
  the L2 ball of radius `C` - but nothing needs it yet.
- **Symbols are for debug printing only.** 5 to 19 have fixed ones; primes
  beyond get arbitrary distinct pairs in subgroup order, so they are stable only
  within one subgroup. `from_accidentals` refuses symbols that collide with each
  other or with the nominals, sharps, flats and octave digits.
- **An accidental defined as one step** - what ups and downs uses in general - is
  still not available, which is why some equal temperaments get no notation.
- **Accidentals off the derived set: the maths does not care, the names do.**
  `from_accidentals` takes whatever list it is given, and everything downstream -
  `temper`, `spell`, `spellings`, the enharmonics, and all of `NotationOptions` -
  treats the accidentals as an opaque list. Handed Johnston's pair, `81/80` and
  the septimal `36/35`, it builds and spells correctly: `7/4` comes out
  `[2, -2, 1, -1]`, which is `16/9` raised a syntonic comma and lowered a
  `36/35`. What needed doing:
  - **Symbols.** Done: an `Accidental` carries its own symbols.
  - **`keeps_nominals`**, hence the recommendation. Done: nominals are degrees,
    a property of the subgroup, and an accidental's ratio never enters.
  - **The half-apotome bound is unenforced** on a supplied accidental, and
    deliberately so: `27/26` is a common accidental for 13 and breaks it. But
    `spelling_cost` charging seven per mark is only sound because of it, so the
    weights may need revisiting for such notations.
  - Supplying a list: done, as `options_with`, `from_temperament_with` and
    `with_count`.

- **An equal temperament sweep is a coverage test, not a judgement.**
  Temperaments can be arbitrarily bad and most of what a sweep turns up is
  nobody's notation. It catches panics and shows what a rule change moved; for
  whether an answer is sensible, read the named list, and the rank 2 and rank 3
  entries especially.

## Where the code lives

- `notation.rs` - the `Notation` type: `generators`, `temper`, `to_interval`,
  `enharmonics`, `spell`, `spellings`, `note`, `nominal_spellings` and
  `keeps_nominals`, the ranking
  `spelling_cost`, the degree map, and the derivation of a single accidental.
- `notation_options.rs` - `Notation::options` and `with_count`: the subset
  search and its ranking. Tested through `options`.
- `simplify.rs` - `Simplifier`: the comma lattice reduced once, then an exact
  weighted L1 enumeration per tempered interval.
- `temperament.rs`, `primes.rs`, `util.rs` - the temperament, the subgroup, and
  integer vector helpers with no music in them.

Examples: `verify` sweeps the invariants and prints a failure count; `dump`
prints a fingerprint of every notation to be diffed across a change; `notations`
lays out each run; `simplifications` and `spellings` are the two questions;
`spellings_width` sweeps how far `spellings` has to walk;
`simplify` and `miracle` walk one temperament. `cargo bench` times the hot path.
