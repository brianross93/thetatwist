import json
import tempfile
import unittest
from pathlib import Path

from experiments.section12_instantiation import (
    cm_bridge_iw_presentation_solver as producer,
)


class IWPresentationSolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context = producer.load_problem_context()

    def test_six_by_sixteen_job_set_is_exact(self):
        self.assertEqual(set(self.context["chart_bits"]), {f"C{i:04b}" for i in range(16)})
        self.assertEqual(set(self.context["translated_reduced"]), set(range(6)))
        self.assertEqual(
            [len(self.context["translated_reduced"][surface]) for surface in range(6)],
            [162, 108, 162, 108, 108, 162],
        )

    def test_constant_row_reductions_prove_the_same_J(self):
        for surface, record in self.context["generator_reductions"].items():
            self.assertEqual(record["surface"], f"Z_{surface}")
            self.assertEqual(record["original_generator_count"], 162)
            self.assertEqual(
                record["reduced_generator_count"],
                len(self.context["translated_reduced"][surface]),
            )
            self.assertEqual(record["identities"], ["R=U*A", "A=C*R"])
            self.assertTrue(record["reduced_from_original_entries"])
            self.assertTrue(record["original_from_reduced_entries"])

    def test_affine_coordinate_trivializations_are_the_projective_ZXY_basis(self):
        z_chart = producer.factor_linear_monomials(0, 2)
        y_chart = producer.factor_linear_monomials(1, 2)
        zero = (0,) * 8
        u = list(zero)
        v = list(zero)
        u[4] = 1
        v[5] = 1
        self.assertEqual(z_chart, (zero, tuple(u), tuple(v)))
        self.assertEqual(y_chart, (tuple(v), tuple(u), zero))

    def test_first_real_problem_has_the_bound_curvilinear_length(self):
        problem = producer.problem_record(self.context, 0, "C0000")
        self.assertEqual(problem["surface"], "Z_0")
        self.assertEqual(problem["surface_section_labels"], [0, 1])
        self.assertEqual(problem["support_points_in_chart"], 59)
        self.assertEqual(problem["expected_curvilinear_length"], 118)
        self.assertEqual(problem["constant_row_basis_generator_count"], 162)
        self.assertEqual(len(problem["input_polynomials"]["curve_equations"]), 4)
        self.assertEqual(len(problem["input_polynomials"]["surface_equations"]), 2)
        self.assertEqual(
            len(problem["input_polynomials"]["J_row_basis_generators"]), 162
        )
        self.assertEqual(len(problem["input_polynomials_sha256"]), 64)

    def test_singular_job_contains_all_exact_certificates(self):
        problem = producer.problem_record(self.context, 1, "C1111")
        script = producer.singular_script(problem)
        self.assertIn("ideal G=liftstd(J,T_G_from_J,Syz_J,\"std\")", script)
        self.assertIn("division(J,G)", script)
        self.assertIn("resolution SchreyerResolution=sres(G,2)", script)
        self.assertIn(
            f"if (quotient_length!={problem['expected_curvilinear_length']})", script
        )
        self.assertIn("assertZeroMatrix(matrix(G)*matrix(SchreyerFirst)", script)

        equality = producer.singular_equality_script(problem)
        self.assertIn("ideal G=std(J)", equality)
        self.assertIn("division(J,G)", equality)
        self.assertNotIn("sres(G,2)", equality)
        self.assertIn(
            f"if (quotient_length!={problem['expected_curvilinear_length']})",
            equality,
        )

        ambient = producer.singular_ambient_equality_script(problem)
        self.assertIn("ideal H=F+J", ambient)
        self.assertIn("ideal G=std(H)", ambient)
        self.assertIn("division(H,G)", ambient)
        self.assertIn(
            f"if (quotient_length!={problem['expected_curvilinear_length']})",
            ambient,
        )

    def test_chartwise_certifying_subideal_keeps_both_conormal_directions(self):
        problem = producer.certifying_subideal_problem_record(
            self.context, 1, "C0111"
        )
        certificate = problem["certifying_subideal"]
        self.assertEqual(certificate["support_blocks_checked"], 51)
        self.assertEqual(certificate["every_selected_support_block_rank"], 2)
        self.assertEqual(certificate["selected_generator_count"], 2)
        self.assertEqual(
            [row["original_column"] for row in certificate["generator_labels"]],
            [0, 1],
        )
        self.assertEqual(
            len(problem["input_polynomials"]["J0_certifying_generators"]), 2
        )
        ambient = producer.singular_ambient_equality_script(problem)
        self.assertIn("J0_selected_conormal_rank_two_every_support", ambient)
        self.assertIn("equal_length_forces_J0_equals_J_equals_IW", ambient)

    def test_cached_schreyer_script_is_bound_to_the_recorded_basis(self):
        problem = producer.problem_record(self.context, 1, "C1111")
        equality = dict(problem)
        equality["reduced_IW_basis_equals_reduced_J_basis"] = [
            {"index": 1, "polynomial": "u0"},
            {"index": 2, "polynomial": "v0"},
        ]
        script = producer.singular_schreyer_script(equality)
        self.assertIn('attrib(G,"isSB",1)', script)
        self.assertIn("resolution SchreyerResolution=sres(G,2)", script)
        self.assertIn("matrix SchreyerMatrix=matrix(SchreyerFirst)", script)
        self.assertNotIn("matrix(SchreyerFirst)[", script)

    def test_static_manifest_has_all_96_jobs(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            manifest_path, reductions_path = producer.write_static_artifacts(
                self.context, output
            )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["job_count"], 96)
            self.assertEqual(len(manifest["jobs"]), 96)
            self.assertEqual(
                len({(row["surface"], row["chart"]) for row in manifest["jobs"]}),
                96,
            )
            self.assertEqual(manifest["cas_required"]["tested_version"], "4.4.1")
            self.assertTrue(reductions_path.is_file())


if __name__ == "__main__":
    unittest.main()
