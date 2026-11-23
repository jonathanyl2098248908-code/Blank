import json
import subprocess
import unittest
from pathlib import Path


ASSET = Path(__file__).parent / "assets" / "hyperfine.js"


class HyperfineTests(unittest.TestCase):
    def solve_node(self, input_value):
        source = f"""
const {{solve}} = require({json.dumps(str(ASSET))});
const result = solve({json.dumps(input_value)});
process.stdout.write(JSON.stringify(result));
"""
        result = subprocess.run(
            ["node", "-e", source],
            check=True,
            capture_output=True,
            text=True,
        )
        return json.loads(result.stdout)

    def solve_product_node(self, input_value):
        source = f"""
const {{solveSpinCoupledProduct}} = require({json.dumps(str(ASSET))});
const result = solveSpinCoupledProduct({json.dumps(input_value)});
process.stdout.write(JSON.stringify(result));
"""
        result = subprocess.run(
            ["node", "-e", source],
            check=True,
            capture_output=True,
            text=True,
        )
        return json.loads(result.stdout)

    def base_input(self):
        return {
            "schema": "orbital-hyperfine-v1",
            "electronic_states": [
                {"id": "e1", "J": "1/2", "parity": "+", "k": 1, "energy_hartree": 0.0,
                 "basis_coefficients": [{"configuration": "CSF-1", "coefficient": 1.0}]}
            ],
            "nuclear_states": [
                {"id": "n1", "I": "1/2", "parity": "+", "T": 0, "k": 1, "energy_mev": 1.0,
                 "basis_coefficients": [{"configuration": "SM-1", "coefficient": 1.0}]}
            ],
            "target": {
                "F": "0",
                "electronic_J": "1/2",
                "electronic_parity": "+",
                "electronic_k": 1,
                "nuclear_state_id": "n1",
            },
            "interactions": [
                {
                    "rank": 0,
                    "parity": "+",
                    "scale_hartree": 0.25,
                    "source": "GRASP transition integrals + KSHELL reduced matrix elements",
                    "electronic_reduced": [{"bra": "e1", "ket": "e1", "value": 2**0.5}],
                    "nuclear_reduced": [{"bra": "n1", "ket": "n1", "value": 2**0.5}],
                }
            ],
        }

    def test_scalar_coupling_energy_and_total_expansion(self):
        result = self.solve_node(self.base_input())
        expected = 0.25 + 1 / 27.211386245988
        self.assertAlmostEqual(result["energy_hartree"], expected, places=12)
        self.assertEqual(result["F"], "0")
        self.assertEqual(result["parity"], "+")
        self.assertAlmostEqual(result["target_overlap"], 1.0)
        self.assertAlmostEqual(sum(term["coefficient"] ** 2 for term in result["total_wavefunction_expansion"]), 1.0)
        self.assertEqual({term["electronic_configuration"] for term in result["total_wavefunction_expansion"]}, {"CSF-1"})
        self.assertEqual({term["nuclear_configuration"] for term in result["total_wavefunction_expansion"]}, {"SM-1"})

    def test_off_diagonal_coupling_mixes_basis_states(self):
        data = self.base_input()
        data["electronic_states"].append({
            "id": "e2", "J": "1/2", "parity": "+", "k": 2, "energy_hartree": 0.1,
            "basis_coefficients": [{"configuration": "CSF-2", "coefficient": 1.0}],
        })
        data["interactions"][0]["electronic_reduced"].extend([
            {"bra": "e2", "ket": "e2", "value": 2**0.5},
            {"bra": "e1", "ket": "e2", "value": 0.2},
        ])
        result = self.solve_node(data)
        self.assertEqual(result["dimension"], 2)
        probabilities = {item["electronic_state"]:item["probability"] for item in result["components"]}
        self.assertGreater(probabilities["e2"], 0)
        self.assertAlmostEqual(sum(item["probability"] for item in result["components"]), 1)

    def test_spin_one_half_dipole_coupling_matches_analytic_f_splitting(self):
        data = self.base_input()
        data["nuclear_states"][0]["energy_mev"] = 0
        data["interactions"] = [{
            "rank": 1,
            "parity": "+",
            "scale_hartree": 1,
            "source": "Analytic J dot I angular-momentum check",
            "electronic_reduced": [{"bra": "e1", "ket": "e1", "value": (1.5) ** 0.5}],
            "nuclear_reduced": [{"bra": "n1", "ket": "n1", "value": (1.5) ** 0.5}],
        }]
        f_zero = self.solve_node(data)
        data["target"]["F"] = "1"
        f_one = self.solve_node(data)
        self.assertAlmostEqual(f_zero["energy_hartree"], -0.75, places=12)
        self.assertAlmostEqual(f_one["energy_hartree"], 0.25, places=12)

    def test_incomplete_transition_table_is_rejected(self):
        data = self.base_input()
        data["electronic_states"].append({
            "id": "e2", "J": "1/2", "parity": "+", "k": 2, "energy_hartree": 0.1,
            "basis_coefficients": [{"configuration": "CSF-2", "coefficient": 1.0}],
        })
        with self.assertRaises(subprocess.CalledProcessError) as error:
            self.solve_node(data)
        self.assertIn("incomplete", error.exception.stderr)

    def test_non_normalized_source_wavefunction_is_rejected(self):
        data = self.base_input()
        data["nuclear_states"][0]["basis_coefficients"][0]["coefficient"] = 0.5
        with self.assertRaises(subprocess.CalledProcessError) as error:
            self.solve_node(data)
        self.assertIn("not normalized", error.exception.stderr)

    def test_ensdf_spin_coupled_product_needs_no_matrix_data(self):
        result = self.solve_product_node({
            "electronic_state": {"id": "GRASP:J1/2-:k1", "J": "1/2", "parity": "-", "k": 1},
            "nuclear_level": {
                "id": "13C:idx0", "spin": "1/2", "parity": "-",
                "spin_parity_assignment_exact": True, "energy_mev": 0.0,
            },
            "F": "0",
        })
        self.assertEqual(result["parity"], "+")
        self.assertEqual(result["energy_offset_mev"], 0)
        self.assertFalse(result["nuclear_basis_resolved"])
        self.assertAlmostEqual(sum(term["probability"] for term in result["total_wavefunction_expansion"]), 1)
        self.assertEqual(len(result["total_wavefunction_expansion"]), 2)

    def test_ensdf_spin_product_rejects_ambiguous_spin_and_invalid_f(self):
        payload = {
            "electronic_state": {"id": "e", "J": "1/2", "parity": "+", "k": 1},
            "nuclear_level": {
                "id": "n", "spin": "1/2", "parity": "+",
                "spin_parity_assignment_exact": False, "energy_mev": 0.0,
            },
            "F": "0",
        }
        source = f"""
const {{solveSpinCoupledProduct}} = require({json.dumps(str(ASSET))});
solveSpinCoupledProduct({json.dumps(payload)});
"""
        with self.assertRaises(subprocess.CalledProcessError) as error:
            subprocess.run(["node", "-e", source], check=True, capture_output=True, text=True)
        self.assertIn("no exact spin/parity", error.exception.stderr)

        payload["nuclear_level"]["spin_parity_assignment_exact"] = True
        payload["F"] = "2"
        source = f"""
const {{solveSpinCoupledProduct}} = require({json.dumps(str(ASSET))});
solveSpinCoupledProduct({json.dumps(payload)});
"""
        with self.assertRaises(subprocess.CalledProcessError) as error:
            subprocess.run(["node", "-e", source], check=True, capture_output=True, text=True)
        self.assertIn("violates the angular-momentum", error.exception.stderr)


if __name__ == "__main__":
    unittest.main()
