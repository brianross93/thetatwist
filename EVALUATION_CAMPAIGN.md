# Corrected Section 12 evaluation campaign

Source: Eyal Markman, [*Secant sheaves and Weil classes on abelian
varieties*](https://arxiv.org/abs/2509.23403v2), Section 12 and Question 11.4.

## Result

The open problem is smaller than the full 16-dimensional Hodge kernel.

Use only the corrected object

$$
E=
\left[\mathcal O_X\longrightarrow\nu_*\mathcal O_{\widetilde Z}\right]
\otimes\mathcal O_X(3\Theta).
$$

Its Chern character is $\alpha+3\beta$. The older $\alpha+\beta$ object is not
part of this campaign.

The Hodge kernel has two parts.

1. Ten classes are ordinary polarized deformation directions:

   $$
   B_{ii},\qquad B_{ij}+B_{ji}\quad(i<j).
   $$

   These classes have zero obstruction. To see this, deform the polarized
   abelian fourfold. Extend each translation point as a section of the
   abelian scheme. The polarization has one relative theta section. Therefore,
   each translated theta divisor extends. The transverse isolated nodes form
   an etale family after a local base change. The selected nodes and the
   partial normalization also extend. The corrected complex $E$ then extends.

2. Six classes are mixed Poisson--gerbe directions:

   $$
   C_{ij}-3A_{ij}\quad(i<j).
   $$

   These six values are the open calculation. They are not ordinary
   deformations of the polarized theta arrangement.

Thus, the next proof target is exact:

$$
\operatorname{ev}_E(C_{ij}-3A_{ij})=0
\quad\text{for all }1\leq i<j\leq4.
$$

Do not replace these six checks with one check unless the fixed divisor and
node configuration has a chain-level symmetry action. A cohomology action is
not sufficient.

## Finite derived presentation

Set $Z_i=D_i\cap D_{i+1}$, with indices modulo $n=(d+9)/2$. The generic
incidence has the following Cech resolution:

$$
0\longrightarrow\mathcal O_Z\longrightarrow C^0
\longrightarrow C^1\longrightarrow C^2\longrightarrow0.
$$

The terms are:

- $C^0$: one codimension-two Koszul complex for each $Z_i$;
- $C^1$: one codimension-three term for each adjacent component pair, and one
  codimension-four term for each nonadjacent pair;
- $C^2$: one codimension-four correction for each consecutive component
  triple.

For $d=3$, the ledger is:

| Item | Count |
|---|---:|
| Surfaces $Z_i$ | 6 |
| Adjacent pair curves | 6 |
| Pairwise point sets | 9 |
| Triple point corrections | 6 |
| Isolated double-point sets | 3 |
| Isolated nodes | 72 |
| Nodes in the partial normalization | 60 |
| Isolated nodes left glued | 12 |

### Parameter--incidence closure check

The corrected parameters close against the incidence resolution. Call this the
**parameter--incidence closure check**. The symbolic incidence schema passes
this check.

For $d=3$, the printed component parameter gives

$$
n=\frac{d+9}{2}=6,
$$

and the resolution has exactly six surfaces $Z_i$. The printed normalization
parameter gives

$$
r=(d+9)(2d-1)=60,
$$

and the resolution normalizes exactly 60 of its 72 isolated nodes. The
remaining node scheme has length

$$
72-60=12=(d+9)(d-2).
$$

More precisely, the schema contains six surfaces $Z_i$, six adjacent triple-
theta curves $T_i$, six consecutive four-theta schemes $Q_i$ of length 24,
and three isolated four-theta schemes of length 24. The last three schemes
contain the 72 isolated nodes. The selected scheme $R_{\mathrm{keep}}$ has
length 12, and its length-60 complement is normalized.

The Cech--Koszul character of this same resolution is the corrected class
$\alpha+3\beta$. Thus, three layers agree:

1. the exact theta-twist character calculation;
2. the printed values of $n$ and $r$;
3. the explicit component and node incidence resolution.

This is a third, structurally separate realization check. It does not treat
the printed parameters as independent premises. It shows that those parameters
are realized by the exact incidence schema, component by component and stratum
by stratum. The later finite-field run must still instantiate this schema with
explicit divisors and nodes. That run will convert the symbolic census into a
fiber-level geometric certificate.

The partial normalization fits in

$$
0\longrightarrow\mathcal O_Z
\longrightarrow\nu_*\mathcal O_{\widetilde Z}
\longrightarrow Q_R\longrightarrow0,
\qquad \operatorname{length}(Q_R)=60.
$$

The extension maps in this sequence are essential. The Chern character does
not determine them.

After the node choice is fixed, the partial normalization has the smaller
incidence resolution

$$
0\longrightarrow\nu_*\mathcal O_{\widetilde Z}
\longrightarrow\bigoplus_i\mathcal O_{Z_i}
\longrightarrow
\left(\bigoplus_i\mathcal O_{T_i}\right)
\oplus\left(\bigoplus_i\mathcal O_{Q_i}\right)
\oplus\mathcal O_{R_{\mathrm{keep}}}
\longrightarrow\bigoplus_i\mathcal O_{Q_i}
\longrightarrow0.
$$

Here, $Q_i$ is a length-24 intersection of four theta divisors. It is the
triple intersection of three consecutive surface components. The scheme
$R_{\mathrm{keep}}$ is the selected reduced length-12 part of the 72 isolated
nodes that stays glued. The complementary 60 nodes are normalized. Section 12
does not select these nodes. A certificate must record the choice.

## Next computation

Do not construct the abelian fourfold by global elimination in the
$3\Theta$ embedding. Since $h^0(X,3\Theta)=3^4=81$, that model lies in
$\mathbf P^{80}$ and is not the useful computational presentation.

Use finite theta-coordinate linear algebra instead. Resolve each
theta-intersection by its Koszul complex. Represent section multiplication and
restriction maps in compatible theta bases. Represent the selected residual
node scheme $R_{\mathrm{keep}}$ by exact point-evaluation and normalization
maps. The universal theta identities do not determine this selected-node data,
so it must remain an explicit input.

The level-3 theta module has dimension 81. However, an exact order-six
translation is not automatically a level-3 clock action. The theta group of
$3\Theta$ covers $X[3]$, while a point of exact order six is not in $X[3]$.
Use either a compatible level-6 theta structure or the induced space

$$
\bigoplus_{i=0}^{5}H^0\!\left(t_{ig}^*3\Theta\right),
$$

where translation permutes the six summands. Neither choice automatically
linearizes the single corrected object $E$. Do not identify the cyclic indexing
with a level-3 Heisenberg action without this additional construction.

A prime $p\equiv1\pmod{12}$ is useful because the base field then contains the
needed small roots of unity. The congruence alone does not make the required
torsion or theta structure rational. For each selected fiber, verify good
reduction, the required rational theta data, and

$$
6g=0,\qquad 2g\ne0,\qquad 3g\ne0.
$$

Also verify the 72-node incidence ledger and the Frobenius stability of
$R_{\mathrm{keep}}$. If the computation will use cyclic symmetry on the node
data, choose $R_{\mathrm{keep}}$ as two full free six-element orbits and record
them. This is an optional symmetry-preserving choice, not a new requirement for
the object. For a characteristic-zero conclusion, reduce one fixed integral
model at independent good primes. Reconstruct from two primes and verify the
result at a third prime. Unrelated finite-field examples do not support
rational reconstruction.

Build the locally free total complex once. Cache its differentials and compute
the principal-parts Atiyah cocycle once. Contract that common cocycle with all
six mixed classes. Then reduce all six results in one exact basis of
$\operatorname{Ext}^2(E,E)$. Macaulay2 can manage the supplied matrices,
complexes, Hom complexes, and syzygy calculations. It need not eliminate the
global $\mathbf P^{80}$ model.

The cyclic translations can reduce the component and node assembly. They act
trivially on the cohomological $HT^2$ source and therefore do not, by
themselves, merge the six mixed witnesses. Inversion also fixes these
two-factor mixed classes. Reduce the six evaluations only if the chosen
configuration has a proved automorphism with a nontrivial linear action on
$H^1(X)$ and a compatible chain-level action on the total complex. Compute its
orbits on $\{C_{ij}-3A_{ij}\}$ before using such a reduction.

The output must include:

- the divisor and selected-node input data;
- the total differential and its basis order;
- the six exact Ext classes;
- the matrix rank and nullspace;
- hashes for all source and result files.

### Declared outcomes

The campaign records both outcomes before the six evaluations run.

- If all six mixed classes vanish, then the known rank-12 contraction and the
  ten geometric vanishings give
  $\ker(\operatorname{ev}_E)=\ker((- )\mathbin{\lrcorner}\operatorname{ch}(E))$
  for this corrected object. This proves Markman's weaker criterion for this
  object. It does not assert surjectivity.
- If one mixed class is nonzero, then this corrected object with this recorded
  divisor and node configuration fails the weaker criterion. A dual cocycle
  that pairs nontrivially with the class is the failure certificate. This result
  does not cover a different choice of residual nodes unless the computation
  also proves independence from that choice.

An input-integrity failure is not a third mathematical outcome. It means that
the run did not instantiate the declared object and has no verdict. An
all-zero finite-field result proves the statement over that field. A complex
result requires the fixed-model lift or reconstruction described above.

Both branches are complete results for the instantiated object. The software
must report the branch that the exact calculation produces.

This calculation is suitable for Macaulay2. Its official macOS distribution
uses Homebrew. Installation is an implementation task, not a mathematical
assumption.

## Claim boundary

The incidence ledger and the $10+6$ reduction are exact. The six mixed Ext
values are not yet computed. This document does not claim a new case of the
Hodge conjecture.

## Implementation references

- D. Lubicz and D. Robert, [*Computing isogenies between abelian
  varieties*](https://arxiv.org/abs/1001.2016), for exact algebraic theta
  coordinates and theta-group computations;
- the official Macaulay2 documentation for
  [matrices](https://macaulay2.com/doc/Macaulay2-1.20/share/doc/Macaulay2/Macaulay2Doc/html/_matrices.html),
  [chain complexes](https://macaulay2.com/doc/Macaulay2-1.23/share/doc/Macaulay2/Complexes/html/___Making_spchain_spcomplexes.html),
  and [finite fields](https://macaulay2.com/doc/Macaulay2-1.24.11/share/doc/Macaulay2/Macaulay2Doc/html/_finite_spfields.html).
