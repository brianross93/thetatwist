# Section 12 PC handoff

## Scope

This handoff runs the finite-overlap rank-two \(K\)-generator tree
contraction used at the current Markman main-lane frontier. It does not run
the entire Section 12 workbench and does not claim \(K\)-coherence, the final
cross-RHom/Yoneda operator, a verdict on Question 11.4, or a result about the
Hodge conjecture.

The observed compatible runtime is:

```text
Python 3.14.6
NumPy 2.4.6
SciPy 1.17.1
```

From the repository root, create an environment and install the exact Python
packages:

```text
python3.14 -m venv .venv-section12
.venv-section12/bin/python -m pip install -r requirements-section12.txt
```

On Windows, the corresponding interpreter is
`.venv-section12\Scripts\python.exe`.

The active path uses only the Python standard library, NumPy, and SciPy. It
does not invoke Sage, Singular, FLINT, NGSolve, or an external service.

## Run

Inspect a generator-specific plan without launching the large calculation:

```text
.venv-section12/bin/python \
  experiments/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_driver.py \
  --generator e3 --dry-run
```

The transport preflight is optional:

```text
.venv-section12/bin/python \
  experiments/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_driver.py \
  --generator e3 --preflight
```

Run the exact producer:

```text
.venv-section12/bin/python \
  experiments/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_driver.py \
  --generator e3
```

The accepted generator labels are `e0`, `e1`, `e2`, and `e3`. Each label has
its own result directory; a preflight is advisory and is not a prerequisite
for the full producer.

Run the focused handoff regression test with:

```text
.venv-section12/bin/python -m unittest -v \
  tests.test_section12_pc_handoff_manifest
```

## Exact Python closure

The driver import closure consists of these 18 modules:

```text
experiments/section12_instantiation/cm_bridge_global_hs_chart_resolution.py
experiments/section12_instantiation/cm_bridge_hs_all_chart_presentations.py
experiments/section12_instantiation/cm_bridge_hs_border_to_local_frames.py
experiments/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit.py
experiments/section12_instantiation/cm_bridge_hs_extension_from_surface_presentation.py
experiments/section12_instantiation/cm_bridge_hs_graph_koszul_j1_overlap_lifts.py
experiments/section12_instantiation/cm_bridge_hs_lambda_overlap_transport.py
experiments/section12_instantiation/cm_bridge_hs_rank2_k_chain_comparison.py
experiments/section12_instantiation/cm_bridge_hs_rank2_k_generator_overlap_j1.py
experiments/section12_instantiation/cm_bridge_hs_rank2_k_generator_residue_transport.py
experiments/section12_instantiation/cm_bridge_hs_rank2_k_generator_shamash.py
experiments/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction.py
experiments/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_driver.py
experiments/section12_instantiation/cm_bridge_hs_shamash_overlap_lifts.py
experiments/section12_instantiation/cm_bridge_iw_hermite_border_basis.py
experiments/section12_instantiation/cm_bridge_iw_presentation_solver.py
experiments/section12_instantiation/cm_bridge_iw_surface_ring_presentation.py
experiments/section12_instantiation/cm_bridge_six_hs_local_seed_emitter.py
```

Keep `experiments/section12_instantiation/__init__.py`, the matching
`tests/test_section12_<module>.py` files, and `tests/__init__.py` with this
closure.

## Frozen input transfer

A fresh generator run requires these 17 frozen upstream files, totaling
11,212,416 bytes in the audited source workspace. The portable `thetatwist`
handoff includes them at these paths. The parent workbench can ignore the
same result tree, so use the manifest transfer if a checkout does not contain
them:

```text
results/section12_instantiation/cm_bridge_global_hs_chart_resolution/global_hs_chart_resolution_certificate.json
results/section12_instantiation/cm_bridge_global_hs_chart_resolution/degree2_global_generator_coefficients.coo.csv
results/section12_instantiation/cm_bridge_global_hs_chart_resolution/product_theta_atlas.csv
results/section12_instantiation/cm_bridge_global_hs_chart_resolution/product_theta_atlas_edges.csv
results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/six_hs_local_seeds.csv
results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/six_hs_support_K_orbits.csv
results/section12_instantiation/descent_section/descent_section_certificate.json
results/section12_instantiation/independent_incidence_fiber/candidate_section_coefficients.csv
results/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit/K_degree2_action_elements.csv
results/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit/K_stable_product_chart_refinement.csv
results/section12_instantiation/cm_bridge_k_equivariant_edge_audit/k_refined_chart_intersections.csv
results/section12_instantiation/cm_bridge_k_equivariant_edge_audit/k_refined_edge_pullbacks.csv
results/section12_instantiation/cm_bridge_hs_extension_from_surface_presentation/hs_extension_certificate.json
results/section12_instantiation/cm_bridge_hs_extension_from_surface_presentation/Z_2_C0010_hs_extension.npz
results/section12_instantiation/cm_bridge_hs_border_to_local_frames/hs_completed_lci_supports.csv
results/section12_instantiation/cm_bridge_iw_hermite_border_basis/Z_2_C0000.json
results/section12_instantiation/cm_bridge_iw_hermite_border_basis/Z_2_C0010.json
```

Preserve these repository-relative paths when copying them. On the source
PC, emit an input-only transfer manifest:

```text
python3 experiments/section12_instantiation/section12_pc_handoff_manifest.py \
  emit > section12-inputs.manifest
```

Copy the manifest and payload to the destination PC, then verify:

```text
python3 experiments/section12_instantiation/section12_pc_handoff_manifest.py \
  verify section12-inputs.manifest
```

The utility is standard-library-only. It emits sorted
`sha256  bytes  relative/path` records, requires the exact selected inventory,
and rejects symlinks instead of following them.

## Large result transfer

Generator outputs remain separate from Git:

```text
e0  results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction/
e1  results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e1/
e2  results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e2/
e3  results/section12_instantiation/cm_bridge_hs_rank2_k_generator_tree_contraction_e3/
```

Completed generator directories currently contain a roughly
0.73--0.98 GB compressed support-map artifact. Transfer only the generator
directories that are needed, using the same external archive or drive as the
frozen inputs. Include chosen result directories in the manifest by repeating
`--generator` during both emission and verification:

```text
python3 experiments/section12_instantiation/section12_pc_handoff_manifest.py \
  emit --generator e2 --generator e3 > section12-e2-e3.manifest

python3 experiments/section12_instantiation/section12_pc_handoff_manifest.py \
  verify section12-e2-e3.manifest --generator e2 --generator e3
```

Do not use `git add -f` for these result directories.
