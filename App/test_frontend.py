import subprocess
import unittest
from pathlib import Path


class FrontendPhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        html = (Path(__file__).parent / "index.html").read_text(encoding="utf-8")
        cls._html = html
        start = html.index("function electronOccupations")
        end = html.index("function makeGraspStarter", start)
        cls.occupation_source = html[start:end]
        start = html.index("function radialScf")
        end = html.index("function createSolverWorker", start)
        cls.radial_scf_source = html[start:end]
        start = html.index("function parseTwoJ")
        end = html.index("function makeGraspStarter", start)
        cls.j_source = html[start:end]
        start = html.index("function parseGraspOutputs")
        end = html.index("function showResearchError", start)
        cls.grasp_output_source = html[start:end]

    def run_node(self, source):
        result = subprocess.run(
            ["node", "-e", source],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()

    def test_all_neutral_occupations_have_the_correct_electron_count(self):
        node = f"""
"use strict";
{self.occupation_source}
let failures = [];
for (let z = 1; z <= 118; z++) {{
  for (const charge of [0, 1, -1]) {{
    if (z - charge < 1) continue;
    try {{
      const occupied = electronOccupations({{z, charge, excitation:"ground"}});
      const count = occupied.reduce((sum, orbital) => sum + orbital.occ, 0);
      if (count !== z - charge) failures.push(`${{z}}/${{charge}}`);
    }} catch (error) {{
      failures.push(`${{z}}/${{charge}}:${{error.message}}`);
    }}
  }}
}}
process.stdout.write(JSON.stringify(failures));
"""
        failures = __import__("json").loads(self.run_node(node))
        self.assertEqual(failures, [])

    def test_hydrogen_and_helium_preview_conserve_electrons(self):
        node = f"""
"use strict";
{self.occupation_source}
{self.radial_scf_source}
for (const z of [1, 2]) {{
  const result = radialScf({{z, charge:0, mass:z === 1 ? 1 : 4, excitation:"ground"}}, () => {{}});
  const count = result.integratedElectrons;
  if (!Number.isFinite(count) || Math.abs(count - z) > 2e-3) {{
    throw new Error(`Z=${{z}} integrated ${{count}} electrons`);
  }}
}}
process.stdout.write("ok");
"""
        self.assertEqual(self.run_node(node), "ok")

    def test_requested_grasp_state_matches_j_parity_and_level_index(self):
        node = f"""
"use strict";
{self.j_source}
{self.grasp_output_source}
const density = [
  "  1      1/2 -",
  "1.0 0.1 0.1",
  "2.0 0.2 0.2",
  "  2      1/2 -",
  "1.0 0.3 0.3",
  "2.0 0.4 0.4",
].join("\\n");
const parsed = parseGraspOutputs(density, "r P(1s)\\n0.1 1.0\\n0.2 0.5\\n");
const selected = requestedDensityState(parsed.densityStates, {{j:"0.5", parity:"-", k:2}});
if (selected.k !== 2 || selected.j !== "1/2" || selected.parity !== "-" || !selected.label.includes("k=2")) {{
  throw new Error("Requested state was not selected exactly");
}}
let rejected = false;
try {{ requestedDensityState(parsed.densityStates, {{j:"1/2", parity:"+", k:2}}); }}
catch (error) {{ rejected = /does not contain requested/.test(error.message); }}
if (!rejected) throw new Error("Mismatched state was accepted");
process.stdout.write("ok");
"""
        self.assertEqual(self.run_node(node), "ok")

    def test_generated_starter_restricts_j_and_checks_configuration_parity(self):
        starter_start = self._html.index("function makeGraspStarter")
        starter_end = self._html.index("function radialScf", starter_start)
        starter_source = self._html[starter_start:starter_end]
        node = f"""
"use strict";
{self.occupation_source}
{self.j_source}
{starter_source}
const result = makeGraspStarter({{z:1, charge:0, mass:1, excitation:"ground", j:"1/2", parity:"+"}});
if (!result.csfAnswers.includes("\\n1,1\\n")) throw new Error("Starter did not restrict GRASP to 2J=1");
let rejected = false;
try {{ makeGraspStarter({{z:1, charge:0, mass:1, excitation:"ground", j:"1/2", parity:"-"}}); }}
catch (error) {{ rejected = /not the requested/.test(error.message); }}
if (!rejected) throw new Error("Wrong-parity starter was accepted");
process.stdout.write("ok");
"""
        self.assertEqual(self.run_node(node), "ok")


if __name__ == "__main__":
    unittest.main()
