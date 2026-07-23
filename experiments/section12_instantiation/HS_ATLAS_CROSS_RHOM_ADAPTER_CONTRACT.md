# Hartshorne--Serre atlas to cross-RHom adapter contract

## Result

The representative whole-ring frame atlas still does not make the global
cross-Hom matrix numerical. It supplies the source-side frame circuits. The
$K$-chain audit supplies the coordinate substitutions, normalized
\(\mathcal O(4H)\) factor, and the correct comparison equations.

The four standard $K$-generator tree contractions are now complete on their
full finite overlap algebras:

| generator | supports / length | $J_1$ | $J_2/L_1$ | factorized $J_3$ | $L_0$ |
|---|---:|---:|---:|---:|---:|
| $e_0$ | $35/70$ | $491\times687$ | $1966\times2750$ | $6866\times9610$ | $492\times688$ |
| $e_1$ | $31/62$ | $435\times687$ | $1742\times2750$ | $6082\times9610$ | $436\times688$ |
| $e_2$ | $35/70$ | $491\times687$ | $1966\times2750$ | $6866\times9610$ | $492\times688$ |
| $e_3$ | $32/64$ | $449\times687$ | $1798\times2750$ | $6278\times9610$ | $450\times688$ |

Together they contain 133 dual-number support factors and close 1,729
factorized \(D_3/J_3\) block identities
(\(455+403+455+416\)). The separated \(p/q\) counit rows replay on each
complete source row, including the \(q\) columns. One
orientation-coboundary law produces every normalization:

\[
\rho_{K,s}=A_K\frac{u_s}{u_t}J_{8,s},\qquad
q=\frac{\rho}{\det B},\qquad
n=\frac{q}{w_H}.
\]

The $e_0$ normalization is the constant pair \((25,0)\); the other three are
nonconstant overlap units. This is one generator-general law, not four fitted
rules.

The coherence adapter now binds the four generator artifacts to a
deterministic registry of 324 typed $K$ edges, 108 order-three cycles, and
486 commutation squares. It emits exact support identities, tuple-key joins,
and composition templates. It does not yet rebase the rectangular maps into
common bases, multiply numeric edge maps around a cell, certify strict cell
closure, or emit an explicit chain homotopy.

The finite comparisons also do not construct a strict
Laurent/surface-ring lift of the counit rows. That lift remains a separate
open scope between the finite comparisons and the open/formal cross-Hom
endpoint. Declaring the target presentation to be the pullback of the source
and using identity matrices would be base change, not a chosen descent
isomorphism, so this contract does not do that.

## Exact chain comparison

On one oriented overlap, let

\[
A_s=\begin{bmatrix}\lambda_s\\-D_{2,s}\end{bmatrix},\qquad
A_t=\begin{bmatrix}\lambda_t\\-D_{2,t}\end{bmatrix}.
\]

If the source finite-algebra length is \(N_s\), use the stored bases

\[
F_{0,s}=Re_*\oplus R^{7N_s+1},\quad
F_{1,s}=R^{28N_s+6},\quad
F_{2,s}=R^{56N_s+6(7N_s+1)},
\]

and the analogous target bases. The producer must emit

\[
J_1:F_{0,s}/Re_*\to F_{0,t}/Re_*,\quad
J_2:F_{1,s}\to F_{1,t},\quad
J_3:F_{2,s}\to F_{2,t},
\]

together with one row \(a:F_{0,s}/Re_*\to R\). With the already known line
factor \(w\), put

\[
L_0=\begin{pmatrix}w&-a\\0&J_1\end{pmatrix},\qquad L_1=J_2.
\]

The four required identities are

\[
J_1D_{2,s}=D_{2,t}J_2,
\]

\[
\lambda_tJ_2-w\lambda_s=aD_{2,s},
\]

\[
L_0A_s=A_tL_1,
\]

and

\[
J_2D_{3,s}=D_{3,t}J_3.
\]

`whole_overlap_chain_jobs.csv` gives the exact matrix shape for every one of
the 192 product-atlas edges. The ranks vary with the chart; the producer does
not replace them by the representative \(688\times2750\) shape. The existing
supportwise Shamash comparisons are bound as regression data, but are not
relabeled as whole-overlap matrices.

## Frame adapter

For a source frame \(y\) and target frame \(x\), write their cleared
inclusion and projection circuits as

\[
\iota_y=J^{\rm num}_y/e_y,\qquad
\pi_x=P^{\rm num}_x/d_x.
\]

The inter-chart transition is then one matrix product:

\[
G_{(t,x)\leftarrow(s,y)}
=\frac{P^{\rm num}_{t,x}L_0J^{\rm num}_{s,y}}
{d_{t,x}e_{s,y}}.
\]

This is the non-tautological extension of the representative-chart formula
\(P_xJ_y/(d_xe_y)\). Once \(L_0\) and the target frame circuits exist, the
adapter must verify inverse transitions and triple cocycles on every nonempty
overlap. K commutation squares and order-three cycles must close strictly or
carry an explicit chain homotopy.

## From frames to cross-Hom

The frame direction is fixed by

\[
v_t=G_{t\leftarrow s}v_s.
\]

For a map \(X:E\to F\), this gives

\[
X_t=G_FX_sG_E^{-1}.
\]

In the column-major basis
\((E_{00},E_{10},E_{01},E_{11})\), the induced Hom matrix is

\[
G_E^{-T}\otimes G_F.
\]

After the whole-atlas frames are restricted into the already emitted formal
and punctured frames, the final comparison block in degree \(q\) is

\[
\Phi_{EF}(f_U,f_P)
=\gamma_F\,\operatorname{res}_U(f_U)\,\gamma_E^{-1}
-\operatorname{res}_P(f_P).
\]

Use it for \((E,F)=(C,A)\) and \((A,C)\). The existing cross-RHom engine then
forms the derived fibres and their Hom differentials. The dimensions of the
open, formal, and punctured Hom coefficient spaces are outputs of that finite
basis assembly. They are not in the current cache, so the contract leaves
them symbolic instead of assigning fabricated scalar-cover dimensions.

## Current boundary

All four standard $K$ generators have passed their complete rectangular
finite-overlap tests. The invertible Hermite decompositions identify their
overlap algebras with 35, 31, 35, and 32 dual-number factors, respectively,
so these are exact finite-algebra calculations rather than samples. The
generator-general orientation-coboundary law above reconstructs the complete
normalization row in every case.

The adapter has also completed the typed registry boundary. All 133
generator-bound supports are indexed, and all 324 edges, 108 cycles, and 486
squares have deterministic composition templates. A template is not a
coherence verdict: the edge maps still use generator-local rectangular bases.
Common-basis rebasing, actual cell products, strict numeric closure, and any
required chain homotopies remain open.

The finite counit rows have not been promoted to strict
Laurent/surface-ring rows. The representative 219-frame atlas is complete,
but only on `Z_2/C0010`. The formal frames and $q_1$ evaluations are complete,
but the target frame atlases and whole-atlas-to-formal restriction are not.

Accordingly, this tranche emits verified finite generator artifacts and a
typed composition registry, not a numeric cross-RHom adapter. Every missing
coefficient remains unknown. No zero placeholder is introduced.

## Artifacts

Outputs are in
`results/section12_instantiation/cm_bridge_hs_atlas_cross_rhom_adapter_contract/`:

- `whole_overlap_chain_jobs.csv`;
- `cross_rhom_dependency_table.csv`;
- `atlas_cross_rhom_adapter_contract.json`;
- `atlas_cross_rhom_adapter_certificate.json`; and
- `manifest.json`.

The four finite generator artifacts are in:

- `results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction/`;
- `results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e1/`;
- `results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e2/`; and
- `results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e3/`.

The typed registry is in
`results/section12_instantiation/cm_bridge_hs_rank2_k_coherence_adapter/`.

## Reproduce

```bash
python3 \
  experiments/section12_instantiation/cm_bridge_hs_atlas_cross_rhom_adapter_contract.py

python3 \
  experiments/section12_instantiation/cm_bridge_hs_rank2_k_coherence_adapter.py

python3 -m unittest -v \
  tests.test_section12_cm_bridge_hs_atlas_cross_rhom_adapter_contract \
  tests.test_section12_cm_bridge_hs_rank2_k_generator_residue_transport \
  tests.test_section12_cm_bridge_hs_rank2_k_generator_tree_contraction_driver \
  tests.test_section12_cm_bridge_hs_rank2_k_coherence_adapter
```

## Claim boundary

This artifact fixes the next producer, all variable matrix shapes, the four
chain equations, the frame adapter, the Hom conjugation convention, and the
Mayer--Vietoris endpoint. The associated generator artifacts bind exact
localized \(J_1/J_2/L_1\), factorized \(J_3\), \(L_0\), and counit
comparisons for all four standard generators on 133 finite support factors.
The associated adapter binds the exact 324-edge, 108-cycle, and 486-square
registry.

It does not emit a strict Laurent/surface-ring lift, common-basis edge
matrices, numeric cycle or square products, a K-coherence verdict, a target
whole-chart frame atlas, a numeric \(\Phi_{EF}\), a cross-Hom differential,
a Yoneda value, an answer to Question 11.4, or a result on the Hodge
conjecture.
