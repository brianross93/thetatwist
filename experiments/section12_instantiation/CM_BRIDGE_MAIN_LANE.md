# Cohen--Macaulay bridge main lane

## Question

The corrected Section 12 object must retain the character

$$
\alpha+3\beta.
$$

An isolated transverse node gives local bivector rank four. An isolated
quadratic tangency lowers this rank to one, but it does not make the map zero.
This calculation tests the positive-dimensional bridge that was selected as
the main lane.

The result has two parts:

1. the curve bridge removes the last local bivector direction; and
2. the ordinary finite conductor seam that restores the global point ledger
   recreates a nonzero local obstruction.

This separation is important. The first statement is about the local
quasi-isomorphism type. The second statement shows that a correct K-theory
identity is not enough.

## 1. The curve bridge closes the local map

Work in

$$
R=k[a,b,c,d].
$$

Let the two surface branches be

$$
U:(a,b)=0,
\qquad
W:(a,c)=0.
$$

They share the curve

$$
\Gamma:(a,b,c)=0.
$$

Their union ideal is

$$
J=(a,b)\cap(a,c)=(a,bc).
$$

The pair $a,bc$ is a regular sequence. The ideal has the exact resolution

$$
0\longrightarrow R
\xrightarrow{(-bc,a)^t}
R^2\longrightarrow J\longrightarrow0.
$$

At a point of $\Gamma$, the derived fiber is present only in degrees $-1$
and $0$. Therefore,

$$
H^2\operatorname{End}(Li_p^*J)=0.
$$

The first-jet calculation gives the same result. Only the $a$ derivative is
nonzero. Every product of two first-jet operators is zero. Thus, all six
divided commutators vanish.

This is the answer to the local question: the positive-dimensional
Cohen--Macaulay bridge removes the last rank-one detector that survived every
isolated tangency.

## 2. The global character still needs 12 points

For $d=3$, the corrected partial-normalization object has

$$
\operatorname{ch}(E_{60})
=1+3x-\frac32x^2-\frac32x^3+\frac38x^4
=\alpha+3\beta.
$$

If all 72 nodes are normalized, then

$$
\operatorname{ch}(E_{72})
=1+3x-\frac32x^2-\frac32x^3-\frac18x^4.
$$

Since $[\mathrm{pt}]=x^4/24$ on a principally polarized abelian fourfold,

$$
12[\mathrm{pt}]=\frac12x^4.
$$

Hence the exact ledger is

$$
\operatorname{ch}(E_{72})+12[\mathrm{pt}]
=\operatorname{ch}(E_{60}).
$$

The curve bridge solves the local obstruction, but the global object must
still produce this $12[\mathrm{pt}]$ difference.

## 3. The natural finite seam does not work

Put $t=d^m$ and define

$$
N_t=\ker\left(
R/(a,b)\oplus R/(a,c)
\longrightarrow R/(a,b,c,t)
\right),
$$

where the map takes the difference of the two restrictions. This is the
honest intermediate normalization along a divisor of length $m$ on
$\Gamma$. It satisfies

$$
0\longrightarrow N_t
\longrightarrow R/(a,b)\oplus R/(a,c)
\longrightarrow\mathcal O_D\longrightarrow0.
$$

Thus, the object $E_t=[R\to N_t]$ has the required K-class:

$$
[E_t]=[E_{72}]+[\mathcal O_D].
$$

For $m=12$, this is the exact target ledger. However, after unit
cancellation, its minimal local complex has ranks

$$
1\longrightarrow4\longrightarrow5\longrightarrow1.
$$

Its differentials are

$$
d_{-2}=
\begin{pmatrix}
t\\-c\\-b\\a
\end{pmatrix},
$$

$$
d_{-1}=
\begin{pmatrix}
-bc&0&-ct&0\\
a&0&0&-t\\
0&-b&c&0\\
0&a&0&c\\
0&0&a&b
\end{pmatrix},
\qquad
d_0=\begin{pmatrix}0&0&a&b&-c\end{pmatrix}.
$$

The producer verifies both identities

$$
d_{-1}d_{-2}=0,
\qquad
d_0d_{-1}=0
$$

over the integral polynomial ring.

For a reduced seam, $t=d$, the six-column local bivector map has rank four.
For a thick seam, $t=d^m$ with $m>1$, it has rank three. Therefore:

- 12 reduced points give 12 copies of the rank-four local model; and
- one fat point of length 12 gives the rank-three local model.

Neither construction has an all-six-zero local map.

This closes the ordinary finite-colength, rank-one conductor seam. Coupling
the seam into the conductor totalization does not fix it.

## 4. The coupled Hecke ledger is exact, but the extension is closed

The degree correction can be written on the six conductor curves without
changing the lower Chern-character terms. Let

$$
i_j:T_j\hookrightarrow X
$$

be the six curves. Choose line bundles $M_j$ and

$$
L_j=M_j(D_j),
\qquad
\sum_j\deg D_j=12.
$$

The symmetric choice is

$$
(\deg D_0,\ldots,\deg D_5)=(2,2,2,2,2,2).
$$

Grothendieck--Riemann--Roch gives

$$
[i_{j*}L_j]-[i_{j*}M_j]
=\deg(D_j)[\mathrm{pt}].
$$

All terms through $\operatorname{ch}_3$ cancel. Only the point term remains.

Set

$$
M=\bigoplus_j i_{j*}M_j,
\qquad
L=\bigoplus_j i_{j*}L_j.
$$

The proposed coupled double Hecke transform first chooses

$$
q:E_{72}\longrightarrow M
$$

and set

$$
H=\operatorname{Fib}(q).
$$

Then choose an inverse-Hecke extension

$$
\xi\in\operatorname{Ext}^1_X(L,H)
$$

and form

$$
H\longrightarrow E_{\mathrm{br}}
\longrightarrow L\xrightarrow{\xi}H[1].
$$

Its K-class is exact:

$$
\begin{aligned}
[E_{\mathrm{br}}]
&=[E_{72}]-[M]+[L]\\
&=[E_{72}]+12[\mathrm{pt}]\\
&=[E_{60}].
\end{aligned}
$$

The K-class is correct, but the required extension does not exist. The local
classifier makes this exact. Put

$$
A=R/(a,b,c)=k[[d]],
\qquad
J=(a,bc),
$$

and use

$$
q:J\longrightarrow A,
\qquad
q(a)=1,
\quad
q(bc)=0.
$$

Its kernel is

$$
H=(a^2,ab,ac,bc).
$$

The sequence

$$
0\longrightarrow H\longrightarrow J\longrightarrow A\longrightarrow0
$$

generates the exact one-dimensional extension module

$$
\operatorname{Ext}^1_R(A,H)\simeq A.
$$

Thus, every scalar extension over the curve DVR is equivalent to one of

$$
1,\qquad d^m\ (m\geq1),\qquad 0.
$$

The result is exhaustive:

- a unit returns the bridge $J$ and has bivector rank zero, but it carries no
  degree;
- a simple zero $d$ has bivector rank six;
- every higher zero $d^m$ has bivector rank three; and
- the split class also has rank three.

Globally,

$$
\mathcal Ext^1_X(i_*L,H)
\simeq i_*\mathcal Hom_T(L,M).
$$

For $L=M(D)$ this becomes

$$
\operatorname{Ext}^1_X(i_*L,H)
\simeq H^0(T,\mathcal O_T(-D)).
$$

When $D$ is effective and nonzero on the proper conductor curve, this group
is zero. The inverse-Hecke class does not globalize.

The same conclusion holds for every ordinary finite-rank matrix extension.
Smith normal form over $k[[d]]$ splits the matrix into scalar blocks
$d^{m_i}$. A positive colength forces at least one nonunit block. Since the
blocks land in separate endomorphism summands, they cannot cancel. Every
positive-colength ordinary matrix has bivector rank at least three.

This closes the rank-one and ordinary vector-bundle double-Hecke lane.

## 5. Hilbert--Burch removes the detector but fails the character ledger

The smallest presentation outside the scalar and Smith categories is

$$
\Phi_s=
\begin{pmatrix}
-bc&b\\
a&c\\
0&s
\end{pmatrix}:R^2\longrightarrow R^3,
$$

where $s=d^m$ cuts the seam. Its signed maximal minors are

$$
as,\qquad bcs,\qquad -(ab+bc^2).
$$

These three elements have no common irreducible factor. A maximal-minor ideal
of a $3\times2$ matrix has height at most two, and the missing common factor
excludes height one. Therefore, the ideal has height exactly two.

The Hilbert--Burch theorem gives the exact resolution

$$
0\longrightarrow R^2
\xrightarrow{\Phi_s}
R^3\longrightarrow I_s\longrightarrow0.
$$

This remains true when $s=0$ at the seam. The derived fiber always has only
two degrees, so

$$
H^2\operatorname{End}(Li_p^*I_s)=0.
$$

All six local bivector classes vanish for a simple seam and for a fat
length-12 seam.

Away from the divisor, $s$ is a unit. The last row and second column cancel,
and the remaining column is

$$
(-bc,a)^t.
$$

Thus, the stabilized object becomes the original CM bridge away from the
seam. This makes it a useful local amplitude-one mechanism.

However, its support audit rejects it as the successor. For $s=d^m$, the
maximal-minor ideal has four minimal primes:

$$
(a,b),\qquad(a,c),\qquad(d,b),\qquad(d,a+c^2).
$$

Their generic multiplicities are

$$
1,\qquad1,\qquad m,\qquad m.
$$

Therefore, its fundamental cycle is

$$
[V(a,b)]+[V(a,c)]
+m[V(d,b)]+m[V(d,a+c^2)].
$$

The scalar seam adds two surfaces. It does not add $m$ points. The lower
Chern character changes before the desired $\operatorname{ch}_4$ correction
is reached.

There are two independent global checks:

1. the literal choice $\alpha=b$, $\beta=c$ is not homogeneous for theta
   bundles; it would require the line bundle of $a$ to be the square of the
   line bundle of $c$; and
2. every ambient line bundle restricts to the triple-theta curve in degree a
   multiple of

   $$
   \Theta^4=24,
   $$

   so it cannot cut the required degree-two divisor.

More generally, an ambient line-pair stabilization has

$$
\operatorname{ch}(Q-P)=e^{qx}(e^{mx}-1).
$$

Its degree-one term is $mx$. The target $12[\mathrm{pt}]$ correction has no
degree-one term. For six regular effective entries, $m_j\geq0$, and
$\sum m_j=0$ forces every $m_j=0$. Then the whole virtual correction is zero.

The correct verdict is:

- the Hilbert--Burch shape proves that a coupled determinantal presentation
  can preserve local amplitude one; but
- this scalar ambient realization has the wrong surface cycle and the wrong
  global Chern character.

It is a design clue, not the missing object.

## 6. The live object is conductor-wide gluing

Let $S=U\cup_TW$ and let $\nu$ be its normalization. Choose branch bundles
$E_U,E_W$, a vector bundle $Q$ on $T$, and an everywhere-surjective map

$$
E_U|_T\oplus E_W|_T\longrightarrow Q.
$$

Set

$$
M=\ker\left(
\nu_*(E_U\oplus E_W)\longrightarrow i_*Q
\right).
$$

This keeps the defect on the full conductor curve. With the required depth
hypotheses, $M$ is maximal Cohen--Macaulay. Compare two conductor bundles
$Q,Q'$ of equal rank. Their curve-supported rank terms cancel, while

$$
\deg Q'-\deg Q=2
$$

changes only the point term by $2[\mathrm{pt}]$. Six curves can therefore
carry the required $12[\mathrm{pt}]$ without using an ambient scalar seam.

For a unital algebra or ideal, the smallest geometric variant considered was
a rank-two square-zero ribbon/Ferrand conductor
$\mathcal O_T\oplus L$ shared by both branches. Sections 7 and 8 record its
classification and the non-transverse refinement.

The standard rank-two local model has now been computed. Use the difference
map

$$
\mathcal O_U^2\oplus\mathcal O_W^2
\longrightarrow\mathcal O_T^2.
$$

Its kernel is $\mathcal O_S^2$. Since $S=V(a,bc)$ is a complete
intersection, the kernel is maximal Cohen--Macaulay and has the direct-sum
Koszul resolution with ranks

$$
2\longrightarrow4\longrightarrow2.
$$

After tensoring this resolution with the closed-point residue field, the
degree-two endomorphism target has dimension four. The exact divided
first-jet screen there has bivector rank zero. The full module-valued Ext
calculation in Section 9 retains the coefficients $c$ and $b$ and has rank
two, so the residue zero is not a local Ext vanishing.

There is also a formal character identity. If $Q,Q'$ have equal rank and

$$
\deg Q'-\deg Q=2,
$$

then

$$
[M_Q]-[M_{Q'}]=[i_*Q']-[i_*Q].
$$

Equal rank cancels the curve term and the normal Todd correction. The degree
difference contributes $2[\mathrm{pt}]$ per curve and
$12[\mathrm{pt}]$ over six curves.

These two checks have not yet been realized by one object. For the standard
difference kernel, both rank-two branch maps are isomorphisms to $Q$. This
fixes the degree of $Q$ from the branch restrictions. A degree-changing
comparison therefore needs either a degeneracy locus or a larger common
source. Those local kernels are not automatically $\mathcal O_S^2$, and
their six cups must be classified.

That ordinary common-source classification is now complete. If one branch
map is invertible, Smith normal form splits the other map into scalar channels

$$
(d^m,1)\quad\text{or}\quad(1,d^m).
$$

The exact minimal resolution of each nonunit channel has ranks

$$
2\longrightarrow4\longrightarrow2.
$$

For $m=1$, its local bivector rank is three, with nonzero columns

$$
a\wedge b,\qquad a\wedge c,\qquad a\wedge d.
$$

For every $m>1$, the tangential first jet disappears, but the two normal
columns remain and the rank is two. A larger source adds pure-branch kernel
summands. Their Koszul cups have rank one: $a\wedge b$ on $U$ and
$a\wedge c$ on $W$.

If neither branch map is surjective at the closed point, the combined map can
be surjective only because at least one target direction is supplied by one
branch alone. Its pure-branch normal component survives. The minimal
complementary rank-loss model was computed explicitly and again has rank
three for a simple zero and rank two for a higher-order zero.

Therefore, an ordinary common-source conductor kernel passes the
closed-point residue screen only on the standard invertible-difference locus.
That locus fixes the degree, and its full module-valued Ext map is still
nonzero. Every ordinary degree-changing path already fails the weaker
residue screen.

## 7. A fixed-support square-zero ribbon does not close the gap

The natural square-zero thickening of the original transverse support is

$$
K=(a^2,ab,ac,bc)=(a,b)(a,c).
$$

It keeps the two reduced surface components, but its exact ideal resolution
has ranks

$$
1\longrightarrow4\longrightarrow4.
$$

The local Atiyah-square calculation has rank three, with the surviving
normal columns $a\wedge b$, $a\wedge c$, and $b\wedge c$. The other
fixed-support orbit leaves the transverse ideal $(a,bc)$ unchanged; it has
zero local cup rank but carries no new degree. Thus a square-zero layer on the
unchanged transverse support cannot provide both requirements.

This statement is deliberately local and specific. It closes fixed-support
square-zero thickenings. It does not close every possible non-transverse
contact geometry.

## 8. Quadratic contact passes the residue screen

The smallest reduced positive-dimensional model that passes the
closed-point residue screen is

$$
Z=V\!\left(a,b(b-c^2)\right).
$$

Its two components are $U=V(a,b)$ and $W'=V(a,b-c^2)$. They meet in the
ribbon

$$
T^{(2)}=V(a,b,c^2),
\qquad
\mathcal O_{T^{(2)}}\cong k[[d]][\epsilon]/(\epsilon^2).
$$

The total space is reduced and Cohen--Macaulay, with the two surface
components each occurring with multiplicity one. Its complete-intersection
resolution is

$$
0\longrightarrow R
\xrightarrow{(-f,a)^t}R^2
\xrightarrow{(a,f)}R
\longrightarrow\mathcal O_Z\longrightarrow0,
\qquad f=b(b-c^2).
$$

After restriction to the reduced conductor $T$, the first jet of $f$ is zero.
Only the $a$-jet remains, and the six divided first-jet cups vanish in that
restricted screen. The degree-two target there remains one-dimensional, so
the screen is nonvacuous. It is not a full module-valued Ext calculation and
is not used as a successor certificate.
The family

$$
V\!\left(a,b\big((1-t)c+t(b-c^2)\big)\right)
$$

connects the original transverse union at $t=0$ to this pinched contact at
$t=1$ while preserving the two component cycles.

The normalization character is also exact at the formal pair level. For two
equal-rank ribbon line bundles with degree difference two, the curve term and
normal Todd term cancel, leaving $2[\mathrm{pt}]$ per conductor curve and
$12[\mathrm{pt}]$ over all six curves.

The original theta geometry blocks this particular globalization. The local
replacement $b-c^2$ compares a section of class $\Theta$ with one of class
$2\Theta$. Homogenizing it as $s b-c^2$ changes the second surface from
class $\Theta^2$ to class $2\Theta^2$. Alternatively, the coefficient needed
to keep the class $\Theta$ would lie in a line bundle of class $-\Theta$ and
has no nonzero global section. Higher contact orders have the same anti-ample
coefficient problem. The pinched model is therefore a valid local design law,
but not a successor inside the unchanged cyclic theta component classes.

## 9. The parity-doubled conductor passes the residue screen, not Ext

There is one explicit local perfect-complex ansatz that combines the zero
local map with the exact virtual point correction. Let

$$
B=R/(a,b),\qquad C=R/(a,c),\qquad A=R/(a,b,c),
\qquad S=R/(a,bc).
$$

The derived branch intersection is

$$
B\otimes_R^{\mathbf L}C\simeq A\oplus A\epsilon[-1],
\qquad \epsilon^2=0.
$$

Choose rank-two conductor bundles $Q,Q'$ with
$\deg Q'-\deg Q=2$ and put them in opposite cohomological parity:

$$
\mathcal Q^\bullet=Q\oplus Q'[1].
$$

The local homotopy fiber is $S^2\oplus S^2[1]$. A bounded free model has
ranks

$$
2\longrightarrow6\longrightarrow6\longrightarrow2.
$$

Its two differential compositions are exactly zero. After tensoring the
minimal free model with the residue field at the closed point, the
degree-two endomorphism target has dimension 24 and all six divided
first-jet commutators vanish. This is an exact residue-field screen. It is
not a calculation of the full local Ext classes, because the tensor product
kills nilpotent coefficients.

The full calculation corrects the earlier conclusion. For
$S=R/(a,bc)$, the Koszul resolution gives

$$
\operatorname{Ext}_R^*(S,S)
\cong S\otimes\Lambda^*(e_1,e_2),
$$

with zero Hom differential. The contracted Atiyah operators are

$$
T_a=e_1,\qquad T_b=c e_2,\qquad T_c=b e_2,\qquad T_d=0.
$$

Hence the two nonzero bivector columns are

$$
a\wedge b\longmapsto c e_{12},\qquad
a\wedge c\longmapsto b e_{12}.
$$

They have module-valued rank two. The diagonal copies in
$S^2\oplus S^2[1]$ therefore survive as actual Ext classes. This correction
does not affect the genuine bridge ideal $J=(a,bc)$: $J$ has projective
dimension one and no local degree-two endomorphism target, whereas $S=R/J$
has projective dimension two and is a different object.

The normalization source cancels in $K$-theory because of the opposite
parities, and

$$
[\operatorname{Fib}]=i_*([Q']-[Q])=2[\mathrm{pt}]
$$

per conductor curve. Six curves give the required $12[\mathrm{pt}]$ exactly
at the virtual-character level. This ledger does not cancel the two
module-valued local cups.
Upper-triangular derived transition functions

$$
\begin{pmatrix}
g_Q&\epsilon h\\
0&g_{Q'}
\end{pmatrix}
$$

can carry the degree difference without a determinant-zero locus because
$\epsilon^2=0$. This observation concerns the conductor transition itself.
It does not yet supply its attachment to the two fixed surface branches.

That attachment has a sharp compatibility condition. For
$\operatorname{Fib}(N\to i_*Q)$ to have the standard local bridge module,
the branch-to-conductor maps are invertible along the conductor and identify
$Q$ with the branch restrictions. Repeating this with the same fixed branch
source for $Q'$ forces

$$
\deg Q'=\deg Q.
$$

The upper-triangular $\epsilon$ term does not change these diagonal
determinant classes. Separate branch sources could avoid the equality only if
$Q$ and $Q'$ both extend to both surfaces; degree difference two is then an
unproved restriction-image condition. Truncating the parity construction to
one ordinary sheaf restores an equal-rank determinant and the Smith-block
obstruction.

The global descent audit makes this split exact. Abstractly, the six curve
bundles exist on the generic complex conductor network: take
$Q_i=\mathcal O_{T_i}^2$ and
$Q'_i=L_i\oplus\mathcal O_{T_i}$ with
$L_i\in\operatorname{Pic}^2(T_i)$, and choose identity trivializations on
the six finite junction schemes. This is valid conductor-nerve data. It is
not an ambient object on the surfaces. For the displayed same-source standard
fiber, the branch restriction argument above gives a degree contradiction.
A nilpotent off-diagonal transition cannot alter it.

The cached theta adapter also fixes the computational consequence. The old
$e_R=0$ base has 25 terms and 48 theta rows. Six local derived templates have
aggregate ranks $12\to36\to36\to12$. They cannot be appended to the old RHom
matrices without the missing branch attachments: the new cross-Hom blocks
require the object differential, principal parts, RHom, and witness columns
to be rebuilt.

This is a virtual-character ansatz and a residue-field screen, not a local or
global successor object.
Buchweitz--Flenner define semiregularity maps for perfect complexes, and
Markman's displayed Section 12 candidate is itself a two-term object on a
smooth variety. Thus the Atiyah, Hochschild-evaluation, and semiregularity
maps are available in the correct category. What is not yet available is the
all-order deformation theorem under only partial semiregularity. It would
follow from the additional statement that every higher lifting obstruction
lies in the image of the Hochschild evaluation map. The first-order ambient
obstruction has this form; persistence through all Artin thickenings is the
missing lemma.

Even after that lemma, a residue-field zero is insufficient. A replacement
must first make the six classes vanish in the full local endomorphism Ext
module. The global object must then have the six corrected evaluations equal
to zero, retain the ten geometric kernel directions, and have those sixteen
directions equal the full annihilator of its Chern character. The complete
category statement is in `DERIVED_CONDUCTOR_CATEGORY_BRIDGE.md`.

## 10. Module-valued conservation law and the remaining lane

The first coupled Smith candidate is the cone of $d^2\operatorname{id}_S$.
It is the Koszul complete intersection

$$
M=R/(a,bc,d^2).
$$

Its closed-point residue screen again has rank zero, but its actual Ext
classes, in the order $ab,ac,ad,bc,bd,cd$, are

$$
c e_{12},\quad b e_{12},\quad 2d e_{13},\quad 0,
\quad 2cd e_{23},\quad 2bd e_{23}.
$$

They have module-valued rank five. The support of this cone is also
one-dimensional, so it is not a point-only correction.

The direct length-12 Artinian complete intersection

$$
P=R/(a^2,b^2,c,d^3)
$$

has the correct scheme length, but its residue-field rank zero hides six
independent module-valued classes:

$$
4ab e_{12},\quad 2a e_{13},\quad 6ad^2 e_{14},
\quad 2b e_{23},\quad 6bd^2 e_{24},\quad 3d^2 e_{34}.
$$

This example is governed by a general conservation law. A nonzero
finite-length perfect block on a smooth fourfold cannot have all six raw
contractions of $\operatorname{At}^2$ vanish. If they did, then
$\operatorname{At}^2=0$, hence $\operatorname{At}^4=0$, and

$$
\operatorname{ch}_4(P)
=\frac{1}{4!}\operatorname{Tr}(\operatorname{At}(P)^4)=0,
$$

contrary to its nonzero virtual length. The stated object type does not forbid
a zero-dimensional coherent summand. The obstruction is intrinsic. A split
sum $E_{72}\oplus P$ still fails because its Atiyah class is block diagonal;
cross-Ext groups do not create off-diagonal evaluation terms for a split
object.

There is also a sharp chain-level restriction. For a one-way coupled
differential

$$
D=\begin{pmatrix}d_A&h\\0&d_B\end{pmatrix},
$$

the $A$-diagonal block of every first-jet commutator remains
$[\partial_i d_A,\partial_j d_A]$. Such an $h$ can help only if the old class
becomes a boundary in the enlarged endomorphism complex. It does not in the
$d^2$ cone, whose Koszul Hom differential against $M$ is zero.

These results close the split parity model, the $d^2$ Smith cone, and a direct
length-12 point correction. They do not close a genuinely non-Smith,
non-split attachment. Such an object must use new attachment maps whose
off-diagonal Atiyah products alter the diagonal classes while preserving the
point-only K-class. The current minimal global candidate is specified in
`NONSPLIT_CONDUCTOR_EXTENSION_DESIGN.md`; larger bidirectional monads remain
possible if that extension shape fails.

## Status table

| Local object | Defect carrier | Local bivector rank | Result |
|---|---|---:|---|
| transverse node | isolated point | 4 | closed as a failure |
| quadratic tangency | isolated length-four point | 1 | closed as a failure |
| CM curve bridge $(a,bc)$ | full curve | 0 | local success |
| reduced finite seam | point on the bridge | 4 | closed as a failure |
| fat length-12 seam | thick point on the bridge | 3 | closed as a failure |
| scalar Hecke class $d$ | point on the curve | 6 | closed as a failure |
| scalar Hecke class $d^m$, $m>1$ | thick point on the curve | 3 | closed as a failure |
| finite-rank Hecke matrix | curve DVR | at least 3 | closed by Smith form |
| scalar Hilbert--Burch stabilization | ambient divisor plus two extra surfaces | 0 | local zero, rejected character |
| standard rank-two difference kernel $S^2$ | full curve bundle | residue 0; module 2 | closed as a full-Ext failure; degree also fixed |
| equal-rank $Q/Q'$ comparison | full curve bundle | conditional/formal | character ledger pass only |
| ordinary common-source degree change | determinant/extra branch channel | at least 1 | closed as a failure |
| fixed-support square-zero thickening | full conductor | 3, or 0 with no degree carrier | closed in this orbit |
| reduced quadratic-contact ribbon | non-transverse full conductor | residue 0 | residue screen only; original theta globalization rejected |
| split parity-doubled conductor ansatz | derived full conductor | residue 0; module 2 | virtual ledger passes, full Ext and same-source attachment fail |
| coupled $d^2$ Smith cone $R/(a,bc,d^2)$ | one-dimensional nonreduced support | residue 0; module 5 | closed as a full-Ext failure; not point-supported |
| direct length-12 CI $R/(a^2,b^2,c,d^3)$ | isolated nonreduced point | residue 0; module 6 | closed by exact Ext and the finite-length conservation law |
| non-Smith non-split attachment | conductor coupled to point defect | true $\operatorname{Ext}^1$: 1 per point, 12 total | split and nonzero scalar orbits both fail all-six boundary test |
| rank-two surface-syzygy pair | six smooth surface components | local $\mathcal Ext^2=0$ | quotient-valid HS seed and exact $12[\mathrm{pt}]$ ledger pass |
| bidirectionally coupled surface pair | both incident surfaces at a crossing | all six polynomial cups zero | exact local Maurer--Cartan solution; global transgression and residues open |

The main point is precise:

> Every tested ordinary degree-changing conductor sheaf reintroduces a local
> normal cup. The parity-doubled model carries the exact virtual
> $12[\mathrm{pt}]$ correction, but its apparent local zero was only a
> residue-field zero: the full module-valued Ext map has rank two. The coupled
> $d^2$ Smith cone has rank five, and a direct length-12 point block has rank
> six. A split point correction is therefore closed. The remaining derived
> lane must be genuinely non-Smith and non-split, so that attachment terms can
> change the diagonal Atiyah products rather than merely accompany them.

## Corrected cross-extension result

Let $B$ be the cached 25-term theta totalization and $P$ the 16-term Koszul
point block. Restricting both free resolutions first gives a
$121\to90\to36$ derived-fiber complex with $29$ cohomology directions per
point. Those 348 directions are Tor created by base change, not
$\operatorname{Ext}^1_X(B,P)$.

Since $B$ is perfect and $P$ has only $H^2(P)=\mathcal O_R$, the ambient
point-supported calculation is instead

$$
1\longrightarrow6\longrightarrow12,
\qquad
(\operatorname{rank}D_0,\operatorname{rank}D_1)=(1,4).
$$

Thus $\dim\operatorname{Ext}^1_X(B,P)=1$ per retained point and 12 over the
selected point algebra. The natural attachment generates each local line.
The phase-2 producer exhausts its two scalar orbits in the full 41-term End
complex at all 12 points. For the split orbit, the six cups span dimension 6;
for the nonzero orbit, they span dimension 4. Every cup is a nonboundary.

This closes the present one-way $B\to P[1]$ shape as a local repair. It does
not close higher-rank or bidirectional conductor objects, nor a
positive-dimensional Cohen--Macaulay bridge. The exact matrices and ledgers
are in the [corrected Ext result directory](../../results/section12_instantiation/cm_bridge_nonsplit_extension_scan/)
and the [two-orbit phase-2 result directory](../../results/section12_instantiation/cm_bridge_nonsplit_phase2_local_screen/).

## Surface-bundle first-syzygy replacement

The smallest live positive-dimensional replacement is now fixed at the
bundle level. On each smooth surface component $i:S\hookrightarrow X$, choose
rank-two bundles $V_0,V_1$ with equal first Chern class and

$$
c_2(V_1)-c_2(V_0)=2[p]
\quad\text{in }\operatorname{CH}_0(S).
$$

Choose one ambient bundle $F$ surjecting onto both $i_*V_j$ and set
$M_j=\ker(F\to i_*V_j)$. Locally each $M_j$ is two copies of the height-two
ideal $(a,b)$ plus a free summand, so it has projective dimension one and no
local degree-two endomorphism target. GRR gives

$$
\operatorname{ch}(M_1)-\operatorname{ch}(M_0)=2[\mathrm{pt}].
$$

The six surface blocks therefore give the exact $12[\mathrm{pt}]$ correction.
This is not realizable by a virtual resolution split only into theta powers:
vanishing of the first four integer moments forces the fourth moment to be a
multiple of $24$, while the target is $12$. The bundle data must be genuinely
nonsplit.

The bundle pair and the local bidirectional coupling are now both explicit.
At a transverse crossing, one block from each incident surface gives a
$6\to16\to9$ complex for which all 16 ordered first-jet products, and hence
all six divided commutators, vanish over the polynomial ring. This is an
exact local result, not a residue-field screen.

The adjacent curve-boundary restriction layer is now closed. For
$D=T_-\cup T_+\in|2H|$, the boundary-oriented Hartshorne--Serre class
$\delta/(xy^2)$ is locally free at $W$ and restricts to zero on the whole
boundary. Thus $V_1|_D\cong V_0|_D$, and the six component splittings give an
abstract cocycle on the adjacent conductor cycle by conjugating the ambient
$V_0$ transitions. The numerical cycle in the certificate is illustrative,
not a computed geometric transition matrix. Point strata are now typed;
supportwise resolution transitions close from edges through four-cells, and
the relative counit closes on every edge. Higher-cell coherence of the counit
homotopies and finite projective descent remain open.

The independent fibre is not cyclic. Instead of transporting the $Z_2$ seed,
the direct theta calculation now supplies six distinct local schemes $W_s$,
all 486 $K$-orbit lifts, six boundary units, and six completed-local
Hartshorne--Serre cocycles. The local seed stage is complete; the global
theta-chart connecting lift is not.

The global-section part of the Hartshorne--Serre lift is now explicit.  On
the common 16-chart product atlas, every surface has an exact
$162\times1296$ value-and-tangent matrix of rank 161.  Two independent jet
normalizations raise the rank to 163 and produce explicit global sections
$X_s,Q_s$ with the required $(x_s,y_s^2)$ germs.  Their $K$-pullbacks cover
all 81 lifts per surface.  This certifies the local generators at all 486
support lifts; it does not yet certify that their globally generated ideal
has no additional base locus.

The split ambient bundle has also been resolved explicitly by six Koszul
blocks. A common evaluation source is
$\mathcal O_X(-7H)^{\oplus2550}$, giving the contracted local shape
$2\to2554\to2$. The compressed section ring now reaches every required
level from 5 through 10. In particular, the formerly missing level-seven
model is exact: it has all 81 $K$-characters, dimension $7^4=2401$ in each
character sector, and exact six-section multiplication programs for
$6\to7$ and $7\to8$. The surface Koszul presentation

$$
H^0(\mathcal O_X(6H))\longrightarrow
H^0(\mathcal O_X(7H))^{\oplus2}\longrightarrow
H^0(\mathcal O_X(8H))\longrightarrow
H^0(\mathcal O_S(8H))
$$

therefore gives the exact quotient dimension
$4096-2(2401)+1296=590$. This closes the section-ring-level gap; it does not
by itself choose a sheaf-valued open-complement Cech representative.
The companion level-eight quotient producer has emitted all six
$590\times4096$ quotient maps, their right inverses, and all 486 formal
restriction factorizations.  Each quotient has rank 590, every right inverse
satisfies $Q_{8,s}S_{8,s}=I_{590}$, and every formal window satisfies
$B_{8,s,k}=T_{8,s,k}Q_{8,s}$.  The result certificate is present, so these
matrices are completed numeric inputs in this status.

The degree-six non-split source is explicit as well. Every surface has a
300-dimensional $H^0(I_W(6H))$ kernel, for 1,800 exact evaluation columns
total. Their $J_{3,5}$ lift coefficients are emitted at all 486 supports,
and the deterministic division formula gives their all-orders formal lifts.
The nontrivial $q_1$ quotient block satisfies both open-formal frame
identities on all 145,800 section/support instances. The universal inner
$\mathcal O_S(8H)$ Cech operator has shape $15\times13$, rank 13, and an
exact inverse on every actual right-hand side; its six-surface direct sum is
$K$-equivariant. These are exact formal and regular-gauge data. The later
chart calculation emits the exact coefficient/restriction data on 96 charts
and 192 edges. It also proves that the contracted $Q\to U$ sector adds no
rows. The remaining gap is the common multiplicative endpoint and its
Yoneda/residue assembly.

The $Q\to U$ selector has also been normalized in the full outer Cech
complex.  Its 104 directions inject into the internal $U$ Hom complex, but
they lie in a contractible $240$-dimensional $Q\xrightarrow{\pm I}E(Q)$
summand. The surviving quotient has dimension $462-240=222$ and contains
the middle and both cross-Hom channels. The subsequent exact outer
calculation produces the adjacent-overlap map and proves $B_RA=0$. It does
not require an extra selector or $e|_Q$ datum. The live gap is the common
multiplicative endpoint and Yoneda assembly, not another $Q\to U$ choice.

The $K$ coordinate action on the chart cover is also explicit.  The original
16 product charts are not $K$-stable, so the correct comparison uses the
1,024 intersections $U_C\cap k^{-1}U_D$.  All 2,048 translated edge-frame
ratios are emitted as exact rational functions, and the degree-two generator
actions pass 7,776 jet equations.  This also proves that the orbit of the two
normalized generators cannot replace the full ideal: it spans at most 162
directions inside the 1,135-dimensional degree-two jet kernel.

The chart-symbol layer is explicit as well.  The full orbit of the product
atlas has 1,296 distinct charts.  On the original charts and edges, the
translated generators span both conormal directions at every support lift;
Nakayama therefore identifies their ideal with $I_{W_s}$ at all 81 support
stalks per surface.  Exact constant-symbol syzygy bases and restriction maps
are emitted for 96 chart systems and 192 edge systems.  Saturation and full
polynomial-module syzygies remain open away from the support.

Those large syzygies are no longer on the critical path for the extension
class itself.  Positive-twist vanishing gives

$$
\mathrm{Ext}^1_S(I_W,\mathcal O_S(4H))
\cong
\mathrm{Ext}^2_S(\mathcal O_W,\mathcal O_S(4H)),
$$

so the six local orientations uniquely recover the six global
Hartshorne--Serre classes.  Purity also marks the exact stopping point of
that shortcut: the off-diagonal degree-zero maps and the
$H^1(\mathcal Ext^1)$ residues retain data on $U=S\setminus W$.

## Exact cross-RHom engine and stabilization boundary

Let $A=E_{72}$ and let the currently defined correction be the split object

$$
C_{\mathrm{split}}=M_1\oplus M_0[1],
\qquad M_j=\operatorname{fib}(q_j:F\to B_j),
\qquad B_j=i_*V_j.
$$

The two cross-Hom objects do not require explicit bases for the rank-2550
kernels. Exact triangles give

$$
R\!\operatorname{Hom}(M_j,A)
\simeq
\operatorname{Cone}\!\left(
R\!\operatorname{Hom}(B_j,A)
\longrightarrow R\!\operatorname{Hom}(F,A)
\right)
$$

and

$$
R\!\operatorname{Hom}(A,M_j)
\simeq
\operatorname{Fib}\!\left(
R\!\operatorname{Hom}(A,F)
\longrightarrow R\!\operatorname{Hom}(A,B_j)
\right),
$$

with the shifts prescribed by $C_{\mathrm{split}}$. The exact finite-field
engine now implements these constructions for two different object
registries using the single law

$$
d_{\mathrm{Hom}}(f)=d_{\mathrm{target}}f-(-1)^h f d_{\mathrm{source}}.
$$

As a regression, it recovers the cached $A\to A$ dimensions
$736,1686,736$, verifies the 25-term/48-arrow $d_A^2=0$ relation, and passes
independent two-object and mapping-cone square-zero tests. No surface-specific
repair parameter occurs in the engine.

There is an important categorical boundary. If both $q_j$ are stabilized by
a zero summand $G$, then

$$
M_j'\simeq M_j\oplus G,
\qquad
C_{\mathrm{split}}'\simeq
C_{\mathrm{split}}\oplus G\oplus G[1].
$$

The added pair has zero K-class but is not acyclic: its differential is zero.
Thus the common source $F$ cannot be cancelled from the split object merely
because its virtual class cancels. A stabilization-independent coupled object
would require new global data, for example a comparison chain map
$m:M_0\to M_1$ with conductor-compatible homotopies, and would be represented
by a relative cone rather than silently identified with
$C_{\mathrm{split}}$. The RHom engine accepts either object registry, but no
such $m$ is claimed here.

The unresolved layer is therefore the compressed global Yoneda calculation.
The existing 48-channel table records fiber-symbol directions only; it is
not a global section basis.  The non-split Hartshorne--Serre calculation has
2,404 coefficients before conductor matching, with surface dimensions
$401,400,401,400,401,401$.  Evaluation into the old 48 fiber directions has
rank 24.  The exact local residue operator has now been emitted with
shape $288\times576$, rank 288, kernel dimension 288, and a sparse right
inverse.

The six base Hodge chains are also no longer missing.  In the standalone
frozen $E_{72}$ endomorphism complex, every

$$
A_{ij}\operatorname{id}_{E_{72}}
$$

is an explicit 25-term diagonal cocycle in the 1,686-dimensional
degree-two basis.  The six columns are independent modulo boundaries, with
an exact dual detector.  This closes the base identity embedding but not its
transport into the coupled successor.

The actual missing package is the common multiplicative open/formal
comparison.  Its degree-zero Cech lifts $b_{\alpha,ij}$ must satisfy

$$
d_{\mathrm{Hom}}b_{\alpha,ij}=h_{\alpha,j}-h_{\alpha,i}.
$$

Their coboundaries supply the Yoneda numerators.  The regular chain
$Q_4\to J_{1,2}\to P_{1,2}\to\mathbf F_{31}$ is already exact at all 486
supports, but it is not the full Mayer--Vietoris formal endpoint.  The typed
registry leaves $\rho_F$, $\rho_U$, the common comparison, and open
cross-channel multiplication absent without zero placeholders.  The typed
principal-part target is not multiplicative: one-sided Laurent boundaries can
have a surviving two-sided product.  Exact finite complexes
$R_x\oplus R_y\to R_{xy}$ now retain those boundary terms for all four frozen
local word families, and all 1,944 projections recover the previous local
matrices.  This closes the strict local coefficient layer but does not create
the still-absent complete open/formal cross-Hom sources.  The typed
factorization is

$$
T=L_{\mathrm{residue}}N_{\mathrm{Yoneda/Cech}}
\operatorname{EvRes}E_{\mathrm{equalizer}}.
$$

For a coherent projective-dimension-one successor, let $T$ be the
transgression from the global surface
$\mathcal Ext^1$ channels to the protected six-dimensional mixed-symbol
space. The exact target is

$$
TH=-3I_6,
$$

and each of the six $H^1(\mathcal Ext^1)$ residue columns must vanish. The
surface blocks have 48 local channels, so there is no local dimension no-go;
their transition functions decide the actual rank of $T$.

The scalar 16-chart nerve has no first cohomology: its two incidence maps
have ranks 15 and 17, and all fifteen edge cocycles now have explicit
primitives.  The adjacent conductor is also exact: it is a shared-normal
cospan, not a full rank-two comparison.  The six curve rows have rank 92
each, so the curve differential has rank 552 and its kernel has dimension
1,852. All 1,852 kernel columns are now explicit in the 2,404-row source
basis. The twelve $84\times194$ upper restriction matrices yield 1,824
columns, and the compatibility kernels yield 28 diagonal columns. The full
embedding has rank 1,852 and passes the six exact curve-row substitution
checks. The live run is therefore
sheaf-valued and point-sensitive. The point layer is already typed as
$D_{\mathrm{point}}=D_{\mathrm{point,raw}}K_{\mathrm{curve}}$; its incidence,
Koszul labels, signs, and retained-projector compression are certified. The
determined branch-labelled precursor has shape $720\times1852$, rank 642,
and kernel dimension 1,210; it is not yet the point differential.

The next input is one universal surface-extension action, not six independent
repair rules and not 300 fitted quotient lifts per surface. On each
geometry-fixed surface, write

$$
J_s=\operatorname{Fib}(M_s\to P_s),\qquad
e_s=(h_s,\omega_s):J_s\to L_s[1],
$$

with

$$
dh_s=-\omega_s\circ r_s.
$$

The abstract class and its completed-local component are fixed. The missing
numeric datum can be supplied by a strict point-aware global/open cochain
$h_s$. The
intrinsic scalar part starts with the $\mathcal O_{Z_s}(4H)$ Koszul model
$16\to162\to256$, not the level-$5\to6$ relation used by the later
$I_W(6H)$ section probe. The full 32-edge atlas remains a valid construction
route, but it is no longer a list of independent prerequisites: it is one
possible realization of this single derived connecting map.

The smallest two-generator chart compression is now excluded on
`Z_2/C0010`. It has two extra simple rational zeros, so 274 of 297 exact
Hermite border relations are nonmembers. A third existing section removes
the two visible points but is still not an equality certificate.  The full
162-generator translated-section ideal does close: degree 11 contains 293
border generators, and a fixed degree-12 prefix of 31,751 Macaulay rows
contains the remaining four.  Since every translated generator already
vanishes with the required value-and-tangent jets, the reverse inclusion is
independent.  Thus

$$
K=I_W
$$

on `Z_2/C0010`.  This is border containment, not Hermite-kernel saturation,
and it does not assert equality on the other 95 charts.

That exact border ideal now has a finite surface-ring first presentation on
the same chart:

$$
S^{1887}\longrightarrow S^{170}\longrightarrow I_W\longrightarrow0.
$$

It comes from the graph-Koszul resolution, six explicit change-of-rings
cycles, and exact unit Tietze contractions.  It is not a length-one
resolution: the left map has higher syzygies.

The representative-chart Hartshorne--Serre cocycle is now exact.  On
`Z_2/C0010`, a graph-Koszul/Shamash carrier has ranks

$$
(1,687,2750,9610),
$$

and its canonical row satisfies, support by support in the polynomial ring,

$$
\lambda d_3=\sum_{j=0}^{5}(-1)^j f_j\mu_j.
$$

Hence it is closed over $S=P/(f_0,\ldots,f_5)$ and represents
$\operatorname{Ext}^2_S(A_W,S)\cong\operatorname{Ext}^1_S(I_W,S)$.  The
exact factored top row and explicit bottom differential give the chart-local
rank-two module $0\to S\to V'_U\to I_W\to0$.  This is an actual cochain,
not a dimension count, but it is not yet a global bundle.

The full six-surface cubical carrier is also indexed.  It contains 96 charts,
192 overlaps, 486 cells, and 1,296 incidence restrictions.  Every $J_0$ is
emitted in canonical value-and-tangent coordinates, and all 384 stored
standard-basis edge maps satisfy their Hermite factorization.  The higher
resolution edge lifts are now materialized too.  On all 192 original atlas
edges and all 11,665 retained overlap factors, exact Laurent-coefficient
supportwise maps $J_0,J_1,J_2,J_3$ close 81,655 degree-two and 151,645
degree-three component identities.  Their $\Gamma$ and $\Omega$ correction
blocks are explicit sparse data.  This closes pairwise edge compatibility,
not by itself coherent descent.  The square layer is now exact on all 144
surface faces and 8,238 face-support incidences.  Every curve sector is
strict.  The two theta sectors share one axis rule: the three pairs inside
$\{0,1,2\}$ have nonzero defects closed by explicit $K_3$ and $K_4$
homotopies, while the three pairs containing axis 3 are strict.  Cube and
four-cube coherence do not follow from this alone.  The cube layer has now
also closed on all 48 surface cubes and 2,585 cube-support incidences.  The
four curve sectors remain strict.  Both theta sectors need nonzero
$K_4/K_5$ homotopies exactly on the 12 cubes with axes $(0,1,2)$; all 36
cubes containing axis 3 are strict.  The six four-cells retain the mandatory
$A_2^2$ term and are exact on all six surfaces.  All 21,888 typed
$A_2^2$ products vanish, while the
individually nonzero theta terms $A_1A_3$ and $A_3A_1$ cancel.  Consequently
$A_4=0$, and the supportwise resolution carrier is coherent through the full
four-cube atlas in the displayed internal-degree window.

The edgewise $\lambda$ transport is also exact.  Its chart-parity-normalized
duality row transforms with weight $h^4$, as required by
$\mathcal O(4H)$.  One simultaneous $68\times16$ system—not six separate
repairs—solves the $K_2$ and all six equation columns on every retained edge
factor.  Across 192 edges and 11,665 factors, every system has rank 12 and
nullity 4, and all 793,220 scalar residual coordinates vanish.  The stored
free-variables-zero row is a canonical representative of a four-dimensional
affine solution family; the solution itself is not unique. The face
calculation is exact as well: all 144 $\mathcal O(4H)$ line cocycles close,
all 8,238 edge-gauge defects vanish, and the canonical face cochain is zero.
The cube and four-cell equations are then forced zero by carrier degree:
all 15,510 cube and 1,824 four-cell relations close without $\Phi_3$ or
$\Phi_4$. This completes the supportwise cubical counit.

The raw four-term carrier has now been promoted on the representative
`Z_2/C0010` whole ring. Its 49 support frames and 170 split-complement frames
form a 219-open strict rank-two atlas, with all 23,871 pair-transition
circuits governed by one retraction--inclusion formula. Direct replay closes
the 160 active complement frames and all 12,720 active complement inverse
pairs, while the support frames bind actual $d_3$ kernel minors. Thus the
representative chart supplies finite strict dg data for the rank-two module.
The four standard $K$ generators now have independently normalized
rectangular comparisons on their complete finite overlap algebras:

| generator | supports / length | $J_1$ | $J_2/L_1$ | factorized $J_3$ | $L_0$ |
|---|---:|---:|---:|---:|---:|
| $e_0$ | $35/70$ | $491\times687$ | $1966\times2750$ | $6866\times9610$ | $492\times688$ |
| $e_1$ | $31/62$ | $435\times687$ | $1742\times2750$ | $6082\times9610$ | $436\times688$ |
| $e_2$ | $35/70$ | $491\times687$ | $1966\times2750$ | $6866\times9610$ | $492\times688$ |
| $e_3$ | $32/64$ | $449\times687$ | $1798\times2750$ | $6278\times9610$ | $450\times688$ |

The four artifacts contain 133 support factors and close 1,729 factorized
$D_3/J_3$ block identities. The separated $p/q$ counit rows replay on the
complete source rows, including the $q$ columns. One orientation-coboundary
law supplies every normalization:

$$
\rho_{K,s}=A_K\frac{u_s}{u_t}J_{8,s},
\qquad
q=\frac{\rho}{\det B},
\qquad
n=\frac{q}{w_H}.
$$

For $e_0$, $n=(25,0)$ is constant. For $e_1,e_2,e_3$, $n$ is a
nonconstant overlap unit. These are instances of one law rather than four
fitted repair rules, and the Hermite decompositions make each result exact
on its full finite overlap algebra rather than a sample.

The typed coherence adapter now binds these artifacts to 324 edge records,
108 order-three cycle records, and 486 commutation-square records. It emits
support identities and composition templates. It does not yet rebase the
rectangular maps into common bases, multiply the numeric edge maps around a
cell, certify strict cell closure, or emit an explicit chain homotopy.

These calculations do not construct a strict Laurent/surface-ring lift of
the finite counit rows. That lift remains a separate open scope before
whole-atlas frame transport and the numeric cross-RHom endpoint. The next
finite computation is common-basis rebasing followed by numeric cycle and
square composition; it is not another generator-specific rule. See
`HS_RANK2_FULL_CHART_ATLAS.md`,
`HS_RANK2_K_CHAIN_COMPARISON.md`, `HS_RANK2_K_GENERATOR_OVERLAP_J1.md`, and
`HS_RANK2_K_GENERATOR_SHAMASH.md`, together with
the four generator artifact directories and
`results/section12_instantiation/cm_bridge_hs_rank2_k_coherence_adapter/`.

Equivalently, one could compute a simultaneous minimal $A_\infty$ model, but
then the required $m_2$, $m_3$, and higher operations must be retained.
Objectwise Ext groups and cohomology ranks do not determine the nested-cone
Yoneda product.  Contracting $h_s$ without transferring those operations
forgets exactly the multiplicative information the final matrix needs.

Once the strict chart models for $e_s$ are compared coherently on overlaps,
derived global sections and tautological evaluation emit $q_0$ and $q_1$ in
common bases. Functorial pre- and postcomposition then
supplies all four maps used by the two exact triangles above. The cross-RHom
engine can immediately assemble both directions and compute their $H^1$.
Only after that assembly can the Yoneda products, $T$, and the six residue
columns be evaluated. Therefore this tranche is an executable reduction of
the last differential, not a Yoneda or Hodge verdict. See
`MARKMAN_MAIN_LANE_STATUS.md`,
`GLOBAL_HS_CHART_RESOLUTION.md`, `AMBIENT_FIRST_SYZYGY_CONTRACT.md`,
`HIGHER_THETA_MULTIPLICATION.md`,
`O8_SECTION_RING_EXTENSION.md`, `O8_OPEN_FORMAL_CECH_OPERATOR.md`,
`O8_SURFACE_QUOTIENT.md`,
`MAPPING_FIBER_STABILIZATION_AUDIT.md`,
`DERIVED_GENERATOR_CROSS_RHOM_REDUCTION.md`,
`FORMAL_GLUING_CROSS_RHOM_AUDIT.md`,
`POINT_INTERFACE_PATH_DGA.md`,
`TYPED_PATH_ENDPOINT_REGISTRY.md`,
`E72_H2_IDENTITY_FACE_CHAINS.md`,
`IW_MACAULAY_MEMBERSHIP_PROBE.md`,
`IW_SURFACE_RING_FIRST_PRESENTATION.md`,
`DEGREE2_PAIR_RESIDUAL_AUDIT.md`,
`K_EQUIVARIANT_EDGE_AUDIT.md`,
`HS_CHART_LOCAL_SYZYGY_AUDIT.md`, `HS_LOCALIZATION_BYPASS_AUDIT.md`,
`HS_SHAMASH_OVERLAP_LIFTS.md`,
`HS_SHAMASH_SQUARE_COHERENCE.md`,
`HS_SHAMASH_CUBE_COHERENCE.md`,
`HS_SHAMASH_FOUR_CELL_COHERENCE.md`,
`HS_LAMBDA_OVERLAP_TRANSPORT.md`, `HS_LAMBDA_OVERLAP_TRANSPORT_FAST.md`,
`HS_LAMBDA_FACE_COHERENCE.md`,
`HS_LAMBDA_CUBE_COHERENCE.md`, `HS_LAMBDA_FOUR_CELL_COHERENCE.md`,
`PROJECTED_Q_COMPUTABILITY_AUDIT.md`,
`GLOBAL_CHANNEL_COHOMOLOGY.md`, `SHARED_NORMAL_CONDUCTOR_AUDIT.md`,
`NORMAL_SYZYGY_EQUALIZER.md`,
`HS_POINT_STRATUM_CONTRACT.md`,
`GLOBAL_HS_THETA_COCYCLE_CONTRACT.md`,
`HS_LOCALIZED_FRAME_INPUT_AUDIT.md`,
`IW_HERMITE_BORDER_BASIS.md`,
`FORMAL_LOCAL_LCI_EDGE_LIFTS.md`,
`REFINED_HS_FORMAL_COVER.md`,
`CECH_YONEDA_LIFT_ENGINE.md`,
`HS_CONDUCTOR_SYNCHRONIZATION.md`, `SIX_HS_LOCAL_SEEDS.md`,
`YONEDA_TRANSITION_SCHEMA.md`, `HS_YONEDA_RESIDUE_AUDIT.md`,
`BIDIRECTIONAL_MC_SOLVER.md`, and `SURFACE_SYMBOL_CANCELLATION.md`.

The old explicit cubic residual link is not an input to this lane. Its first
equation mixes different $K$-characters, so it is only a pullback calculation
and does not descend to the quotient. The quotient-valid curvilinear seed and
the original Hartshorne--Serre pair remain valid.

If such a block passes, the theta adapter gives the exact next global run:
construct branch extensions and attachment maps, add their $K$-linearizations
and principal parts, rebuild RHom including cross-Homs, and evaluate the six
corrected witnesses and ten geometric directions. In parallel, the
ambient-image obstruction property in
`DERIVED_CONDUCTOR_CATEGORY_BRIDGE.md` remains the theorem needed to turn
kernel equality into the deformation conclusion sought in Question 11.4.

## Reproduce

Run:

```text
python3 experiments/section12_instantiation/cm_bridge_main_lane.py
python3 experiments/section12_instantiation/cm_bridge_hecke_classifier.py
python3 experiments/section12_instantiation/cm_bridge_hilbert_burch.py
python3 experiments/section12_instantiation/cm_bridge_globalization_audit.py
python3 experiments/section12_instantiation/cm_bridge_conductor_wide.py
python3 experiments/section12_instantiation/cm_bridge_common_source_classifier.py
python3 experiments/section12_instantiation/cm_bridge_ribbon_classifier.py
python3 experiments/section12_instantiation/cm_bridge_pinched_ribbon.py
python3 experiments/section12_instantiation/cm_bridge_pinched_globalization_audit.py
python3 experiments/section12_instantiation/cm_bridge_derived_conductor.py
python3 experiments/section12_instantiation/cm_bridge_derived_descent_audit.py
python3 experiments/section12_instantiation/cm_bridge_derived_theta_adapter.py
python3 experiments/section12_instantiation/cm_bridge_coupled_parity_classifier.py
python3 experiments/section12_instantiation/cm_bridge_nonreduced_ci_atyah_audit.py
python3 experiments/section12_instantiation/cm_bridge_nonsplit_extension_scan.py
python3 experiments/section12_instantiation/cm_bridge_nonsplit_phase2_local_screen.py
python3 experiments/section12_instantiation/cm_bridge_surface_syzygy_pair.py
python3 experiments/section12_instantiation/cm_bridge_six_hs_local_seed_emitter.py
python3 experiments/section12_instantiation/cm_bridge_hs_conductor_synchronization.py
python3 experiments/section12_instantiation/cm_bridge_global_hs_chart_resolution.py
python3 experiments/section12_instantiation/cm_bridge_ambient_first_syzygy_contract.py
python3 experiments/section12_instantiation/cm_bridge_higher_theta_multiplication.py
python3 experiments/section12_instantiation/cm_bridge_o8_section_ring_extension.py
python3 experiments/section12_instantiation/cm_bridge_o8_surface_quotient.py
python3 experiments/section12_instantiation/cm_bridge_q1_open_formal_evaluation.py
python3 experiments/section12_instantiation/cm_bridge_o8_open_formal_cech_operator.py
python3 experiments/section12_instantiation/cm_bridge_mapping_fiber_stabilization_audit.py
python3 experiments/section12_instantiation/f31_cross_rhom.py
python3 experiments/section12_instantiation/cm_bridge_k_equivariant_edge_audit.py
python3 experiments/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit.py
python3 experiments/section12_instantiation/cm_bridge_hs_localization_bypass_audit.py
python3 experiments/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction.py
python3 experiments/section12_instantiation/cm_bridge_hs_rank2_k_coherence_adapter.py
python3 experiments/section12_instantiation/cm_bridge_hs_atlas_cross_rhom_adapter_contract.py
python3 -m unittest \
  tests.test_section12_cm_bridge_main_lane \
  tests.test_section12_cm_bridge_hecke_classifier \
  tests.test_section12_cm_bridge_hilbert_burch \
  tests.test_section12_cm_bridge_globalization_audit \
  tests.test_section12_cm_bridge_conductor_wide \
  tests.test_section12_cm_bridge_common_source_classifier \
  tests.test_section12_cm_bridge_ribbon_classifier \
  tests.test_section12_cm_bridge_pinched_ribbon \
  tests.test_section12_cm_bridge_pinched_globalization_audit \
  tests.test_section12_cm_bridge_derived_conductor \
  tests.test_section12_cm_bridge_derived_descent_audit \
  tests.test_section12_cm_bridge_derived_theta_adapter \
  tests.test_section12_cm_bridge_coupled_parity_classifier \
  tests.test_section12_cm_bridge_nonreduced_ci_atyah_audit \
  tests.test_section12_cm_bridge_nonsplit_extension_scan \
  tests.test_section12_cm_bridge_nonsplit_phase2_local_screen \
  tests.test_section12_cm_bridge_surface_syzygy_pair \
  tests.test_section12_cm_bridge_six_hs_local_seed_emitter \
  tests.test_section12_cm_bridge_hs_conductor_synchronization \
  tests.test_section12_cm_bridge_global_hs_chart_resolution \
  tests.test_section12_cm_bridge_ambient_first_syzygy_contract \
  tests.test_section12_cm_bridge_higher_theta_multiplication \
  tests.test_section12_cm_bridge_o8_section_ring_extension \
  tests.test_section12_cm_bridge_o8_surface_quotient \
  tests.test_section12_cm_bridge_q1_open_formal_evaluation \
  tests.test_section12_cm_bridge_o8_open_formal_cech_operator \
  tests.test_section12_cm_bridge_mapping_fiber_stabilization_audit \
  tests.test_section12_f31_cross_rhom \
  tests.test_section12_cm_bridge_k_equivariant_edge_audit \
  tests.test_section12_cm_bridge_hs_chart_local_syzygy_audit \
  tests.test_section12_cm_bridge_hs_localization_bypass_audit \
  tests.test_section12_cm_bridge_hs_rank2_k_generator_tree_contraction \
  tests.test_section12_cm_bridge_hs_rank2_k_generator_residue_transport \
  tests.test_section12_cm_bridge_hs_rank2_k_generator_tree_contraction_driver \
  tests.test_section12_cm_bridge_hs_rank2_k_coherence_adapter \
  tests.test_section12_cm_bridge_hs_atlas_cross_rhom_adapter_contract -v
```

The exact artifacts are in:

- `results/section12_instantiation/cm_bridge_main_lane/`
- `results/section12_instantiation/cm_bridge_hecke_classifier/`
- `results/section12_instantiation/cm_bridge_hilbert_burch/`
- `results/section12_instantiation/cm_bridge_globalization_audit/`
- `results/section12_instantiation/cm_bridge_conductor_wide/`
- `results/section12_instantiation/cm_bridge_common_source_classifier/`
- `results/section12_instantiation/cm_bridge_ribbon_classifier/`
- `results/section12_instantiation/cm_bridge_pinched_ribbon/`
- `results/section12_instantiation/cm_bridge_pinched_globalization_audit/`
- `results/section12_instantiation/cm_bridge_derived_conductor/`
- `results/section12_instantiation/cm_bridge_derived_descent_audit/`
- `results/section12_instantiation/cm_bridge_derived_theta_adapter/`
- `results/section12_instantiation/cm_bridge_coupled_parity_classifier/`
- `results/section12_instantiation/cm_bridge_nonreduced_ci_atyah_audit/`
- [`results/section12_instantiation/cm_bridge_nonsplit_extension_scan/`](../../results/section12_instantiation/cm_bridge_nonsplit_extension_scan/)
- [`results/section12_instantiation/cm_bridge_nonsplit_phase2_local_screen/`](../../results/section12_instantiation/cm_bridge_nonsplit_phase2_local_screen/)
- [`results/section12_instantiation/cm_bridge_surface_syzygy_pair/`](../../results/section12_instantiation/cm_bridge_surface_syzygy_pair/)
- [`results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/`](../../results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/)
- [`results/section12_instantiation/cm_bridge_hs_conductor_synchronization/`](../../results/section12_instantiation/cm_bridge_hs_conductor_synchronization/)
- [`results/section12_instantiation/cm_bridge_global_hs_chart_resolution/`](../../results/section12_instantiation/cm_bridge_global_hs_chart_resolution/)
- [`results/section12_instantiation/cm_bridge_ambient_first_syzygy_contract/`](../../results/section12_instantiation/cm_bridge_ambient_first_syzygy_contract/)
- [`results/section12_instantiation/cm_bridge_higher_theta_multiplication/`](../../results/section12_instantiation/cm_bridge_higher_theta_multiplication/)
- [`results/section12_instantiation/cm_bridge_o8_section_ring_extension/`](../../results/section12_instantiation/cm_bridge_o8_section_ring_extension/)
- [`results/section12_instantiation/cm_bridge_q1_open_formal_evaluation/`](../../results/section12_instantiation/cm_bridge_q1_open_formal_evaluation/)
- [`results/section12_instantiation/cm_bridge_o8_open_formal_cech_operator/`](../../results/section12_instantiation/cm_bridge_o8_open_formal_cech_operator/)
- [`results/section12_instantiation/cm_bridge_mapping_fiber_stabilization_audit/`](../../results/section12_instantiation/cm_bridge_mapping_fiber_stabilization_audit/)
- [`results/section12_instantiation/f31_cross_rhom/`](../../results/section12_instantiation/f31_cross_rhom/)
- [`results/section12_instantiation/cm_bridge_k_equivariant_edge_audit/`](../../results/section12_instantiation/cm_bridge_k_equivariant_edge_audit/)
- [`results/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit/`](../../results/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit/)
- [`results/section12_instantiation/cm_bridge_hs_localization_bypass_audit/`](../../results/section12_instantiation/cm_bridge_hs_localization_bypass_audit/)
- [`results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction/`](../../results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction/)
- [`results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e1/`](../../results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e1/)
- [`results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e2/`](../../results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e2/)
- [`results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e3/`](../../results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e3/)
- [`results/section12_instantiation/cm_bridge_hs_rank2_k_coherence_adapter/`](../../results/section12_instantiation/cm_bridge_hs_rank2_k_coherence_adapter/)
- [`results/section12_instantiation/cm_bridge_hs_atlas_cross_rhom_adapter_contract/`](../../results/section12_instantiation/cm_bridge_hs_atlas_cross_rhom_adapter_contract/)

## References for the category distinction

- [Eyal Markman, *Secant sheaves and Weil classes on abelian varieties*](https://arxiv.org/abs/2509.23403)
- [Ragnar-Olaf Buchweitz and Hubert Flenner, *A semiregularity map for modules and applications to deformations*](https://arxiv.org/abs/math/9912245)
