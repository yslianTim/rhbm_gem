import unittest

from joint_fixed_neighbor_core_size import analyze


def row(topology, atoms, core, search, converged=True, kkt=1e-12, width=1e-14):
    return {"topology": topology, "atoms": atoms, "core_atoms": core,
            "policy": "OneAccepted", "status": "completed", "search_converged": converged,
            "final_global_ac_kkt": kkt, "final_raw_width_gradient_inf_norm": width,
            "search_seconds": search, "sweeps": 7, "confirmed_stationarity_sweep": 7}


class CoreSizeTest(unittest.TestCase):
    def test_selects_fastest_core_when_all_points_are_correct(self):
        rows = [row("chain", 512, 64, 30.), row("chain", 512, 128, 20.), row("chain", 512, 256, 25.),
                row("cube", 512, 64, 40.), row("cube", 512, 128, 35.), row("cube", 512, 256, 45.)]
        rows += [row("cube", 1024, 64, 80.), row("cube", 1024, 128, 70.), row("cube", 1024, 256, 90.)]
        report = analyze(rows)
        self.assertEqual(report["correctness_gate"], "passed")
        self.assertEqual(report["selected_core_size"], 128)

    def test_failed_core_is_not_selected(self):
        rows = [row("chain", 512, 64, 30.), row("chain", 512, 128, 20.),
                row("chain", 512, 256, 25., converged=False)]
        report = analyze(rows)
        self.assertEqual(report["correctness_gate"], "failed")
        self.assertIsNone(report["selected_core_size"])


if __name__ == "__main__":
    unittest.main()
