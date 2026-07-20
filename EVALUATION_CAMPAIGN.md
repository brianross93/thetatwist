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

Build a locally free total complex from the Koszul and Cech terms. Add the 60
point-normalization extension maps. Compute the principal-parts Atiyah cocycle
of the total differential. Contract that cocycle with each of the six mixed
classes. Reduce each result in an exact basis of $\operatorname{Ext}^2(E,E)$.

The output must include:

- the divisor and selected-node input data;
- the total differential and its basis order;
- the six exact Ext classes;
- the matrix rank and nullspace;
- hashes for all source and result files.

This calculation is suitable for Macaulay2. Its official macOS distribution
uses Homebrew. Installation is an implementation task, not a mathematical
assumption.

## Claim boundary

The incidence ledger and the $10+6$ reduction are exact. The six mixed Ext
values are not yet computed. This document does not claim a new case of the
Hodge conjecture.
