# Notation and temperament

Everything below is elementary algebra over `Z`. It is written out because the
statements are useful and easy, not because they are hard.

## Setup

Fix primes `p_1 < ... < p_n`. The **intervals** over them form the free abelian
group `Z^n`, written multiplicatively: the vector `(a_1, ..., a_n)` is the
positive rational `p_1^a_1 * ... * p_n^a_n`. Addition of vectors is
multiplication of rationals.

A **temperament** is a surjective homomorphism

    T : Z^n -> Z^r.

Musically it declares which intervals are to be treated as the same: two
intervals are the same pitch exactly when they have the same image. Its kernel
`ker T` is the lattice of **commas**, the intervals tempered to a unison, and it
has rank `n - r`. (Surjectivity is the usual "no contorsion" condition. Without
it `Z^n / ker T` has torsion, and some interval would be a unison only after
being stacked up several times.)

A **notation** is a surjective homomorphism

    N : Z^n -> Z^m

together with a choice of basis of `Z^m`. Musically the basis elements are the
octave, the fifth, and one accidental apiece for the remaining primes, and a
vector in `Z^m` is a written note. Its kernel `ker N` is the lattice of
intervals the notation writes identically, of rank `n - m`.

The two are homomorphisms of the same shape. The whole subject is what happens
when you ask them to agree.

## 1. When a notation notates a temperament

The following are equivalent.

1. `ker N` is contained in `ker T`.
2. There is a homomorphism `t : Z^m -> Z^r` with `T = t . N`, and it is unique.
3. Any two intervals with the same spelling are the same pitch.

*Proof.* (1) implies (2): `N` is surjective, so it identifies `Z^n / ker N` with
`Z^m`; a homomorphism out of `Z^n` killing `ker N` descends to that quotient,
and the descent is unique because `N` is onto. (2) implies (1) at once, and (3)
is a restatement of (1). Note `t` is surjective, since `T` is. ∎

Condition (1) is the whole compatibility requirement, and the inclusion goes
only one way. A notation is allowed to spell fewer things alike than the
temperament calls alike; it is never allowed to spell more. The gap between the
two kernels is what the next section is about.

## 2. The enharmonic lattice

Define

    E := ker t,

the intervals of notation that the temperament, having factored through, still
sends to nothing. Equivalently `E = N(ker T)`: if `c` is a comma then
`t(N(c)) = T(c) = 0`; and conversely any `y` in `ker t` is `N(x)` for some `x`,
which then has `T(x) = t(y) = 0`.

These are the **enharmonics**: pairs of written notes that the temperament calls
one pitch. `C#` and `Db` in twelve-tone equal temperament differ by an element
of `E`.

`E` is free, being a subgroup of a free abelian group, and it is saturated in
`Z^m`, since `Z^m / E` embeds in `Z^r` and so is torsion-free.

## 3. Commas and enharmonics account for everything

**Theorem.** The sequence

    0 -> ker N -> ker T -> E -> 0

is exact, the second map being inclusion and the third the restriction of `N`.
Hence `ker T / ker N` is isomorphic to `E`, the sequence splits, and

    ker T  ≅  ker N  ⊕  E.

*Proof.* Inclusion is fact 1. `N` carries `ker T` onto `E` by the description of
`E` above. The kernel of `N` restricted to `ker T` is `ker N ∩ ker T = ker N`,
again by fact 1. A short exact sequence whose right-hand term is free splits. ∎

Concretely: take any commas `c_1, ..., c_k` of the temperament whose spellings
`N(c_1), ..., N(c_k)` form a basis of `E`. Together with a basis of `ker N` they
form a basis of `ker T`. Nothing is left over and nothing is counted twice.

**Corollary (ranks).** `rank(ker N) + rank(E) = rank(ker T)`, which is

    (n - m) + (m - r) = n - r.

So the tempered-out dimensions of a temperament are split cleanly in two by any
notation of it: those the notation also spells alike, and those it spells apart.
Every dimension is one or the other.

## 4. When spelling is faithful

The following are equivalent: `E = 0`; `m = r`; `ker N = ker T`; every pitch has
exactly one spelling. In that case `N` induces an isomorphism
`Z^n / ker T ≅ Z^m`, and notation and temperament are the same map in different
clothes.

## 5. An equal temperament always has enharmonics

Whenever `r < m` we have `rank E = m - r > 0`, so `E` is not trivial. A notation
has at least an octave and a fifth among its basis elements, so `m >= 2`, and
this applies to every temperament of rank one.

Take then `r = 1`: an **equal temperament**. Some multiple of the
fifth equals some multiple of the octave there — the circle of fifths closes.
But `Z^m` is free on its basis, and distinct basis elements of a free abelian
group satisfy no relation whatever. So the circle of fifths cannot close in
notation, and the relation that closes it in the temperament is an enharmonic.
In twelve-tone equal temperament that relation is twelve fifths against seven
octaves, and the enharmonic is the Pythagorean comma, which is exactly what
leaves `C#` and `Db` apart on the page.

This is not a defect of any particular notation. Freeness is what makes a
notation writable at all: you want to be able to write down any stack of
generators and have it mean something.

## 6. Three cosets

Fix a pitch, that is, an element of `Z^r` in the image of `T`.

- The just intervals realizing it are a coset of `ker T` in `Z^n`.
- Its spellings are a coset of `E` in `Z^m`.
- The just intervals sharing one spelling are a coset of `ker N` in `Z^n`.

Fact 3 says exactly that the first is fibred over the second by `N`, with fibres
the third. Choosing a canonical member of any of these cosets — the simplest
ratio for a pitch, or the easiest spelling — is a closest-vector problem on a
lattice, once a notion of size is fixed.

Two remarks on that. The natural measures of complexity here are `L1` norms in
the prime-exponent coordinates: the sum of the prime factors of numerator and
denominator, with multiplicity, or the same with each prime replaced by its
logarithm. Lattice algorithms want a quadratic form, so the `L1` problem has to
be solved some other way. And the norm must be expressed in the prime
coordinates, where each coordinate has a prime attached to it; there is no
meaningful way to weight the coordinates of `Z^m`, since a fifth is not a prime.

## 7. Where a prime sits on the fifth chain

A lemma about notations of the kind described above, where each prime beyond the
first two is written with its own accidental.

Suppose the accidental for a prime `p` is

    a = 2^α * 3^β * p^s,    s = ±1,

which is the usual shape: a small interval correcting a point of the fifth chain
to `p`. Give `Z^m` the basis in which the octave is `2` and the fifth is `3/2`,
and write the **fifth coordinate** of an interval for its second coordinate in
that basis — how far along the chain of fifths it lies.

**Lemma.** The fifth coordinate of `p` is `-s * β`.

*Proof.* From `a = 2^α 3^β p^s` and `s² = 1` we get `p = a^s 2^{-sα} 3^{-sβ}`.
The accidental `a` is a basis element other than the octave and the fifth, so
its fifth coordinate is `0`; the octave `2` has fifth coordinate `0`; and
`3 = 2 * (3/2)` has fifth coordinate `1`. Fifth coordinates add, so the fifth
coordinate of `p` is `-sβ`. ∎

For the syntonic comma `81/80 = 2^-4 3^4 5^-1` this gives `s = -1`, `β = 4`, so
`5` sits four fifths up — on `E`, counting from `C`. For `64/63` it puts `7` two
fifths down, and for `33/32` it puts `11` one fifth down.

In Western notation the seven letter names repeat every seven fifths, since
seven fifths is a sharp and a sharp does not change the letter. So the letter a
prime is written on is this number modulo `7`, and two notes share a letter
exactly when their fifth coordinates agree modulo `7`.
