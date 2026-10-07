import unittest

import joint_optimized_frontier as frontier


def record(route, seconds, formal=True):
    result = {"status": "completed" if formal else "timeout", "completed": formal,
              "search_seconds": seconds if formal else None, "total_seconds": seconds,
              "peak_rss_mb": 100., "search_completed": formal}
    return {
        "topology": "chain", "atoms": 512, "rows": 1000, "parameter_count": 1536,
        "route": route, "measurement_scope": "search-only", "input_sha256": "same",
        "formal_result": result, "diagnostic_result": result,
    }


class OptimizedFrontierTest(unittest.TestCase):
    def test_optimized_route_wins_when_it_is_the_fastest_majority(self):
        records = []
        for atoms, optimized_seconds in ((512, 10.), (768, 12.), (1024, 14.)):
            for route, seconds in (("LegacyCompact", 80.), ("OperatorPcg", 40.),
                                   ("FixedNeighbor", 60.),
                                   ("FixedNeighborOptimized", optimized_seconds)):
                row = record(route, seconds)
                row["atoms"] = atoms
                records.append(row)
        analysis = frontier._analyze(records)
        self.assertEqual(analysis["route_position"], "FixedNeighborOptimized-time-dominant")
        self.assertEqual(len(analysis["groups"]), 3)
        self.assertEqual(analysis["groups"][0]["fastest_route"], "FixedNeighborOptimized")

    def test_operator_position_is_preserved_when_optimized_is_slower(self):
        records = []
        for atoms in (512, 768):
            for route, seconds in (("LegacyCompact", 80.), ("OperatorPcg", 10.),
                                   ("FixedNeighbor", 60.),
                                   ("FixedNeighborOptimized", 40.)):
                row = record(route, seconds)
                row["atoms"] = atoms
                records.append(row)
        self.assertEqual(frontier._analyze(records)["route_position"], "OperatorPcg-dominant")

    def test_parser_defaults_to_qualified_policy(self):
        args = frontier.build_parser().parse_args(["--build-dir", "build", "--output-dir", "out"])
        self.assertEqual(args.fixed_local_work, "one")
        self.assertEqual(args.fixed_core_atoms, 64)


if __name__ == "__main__":
    unittest.main()
