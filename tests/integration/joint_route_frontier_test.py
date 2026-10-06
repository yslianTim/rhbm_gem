import csv
import tempfile
from pathlib import Path
from unittest.mock import patch

import joint_route_frontier as frontier


def record(topology, atoms, route, completed, seconds, rss, scope="search-only"):
    result = {"status": "completed" if completed else "timeout", "completed": completed,
              "search_seconds": seconds if completed else None, "total_seconds": seconds,
              "peak_rss_mb": rss, "search_completed": completed}
    return {
        "topology": topology, "atoms": atoms, "rows": 1000 * atoms,
        "parameter_count": 3 * atoms, "route": route,
        "measurement_scope": scope, "formal_result": result,
        "diagnostic_result": result,
    }


def test_frontier_reports_fastest_rss_and_formal_envelope_separately():
    records = [
        record("chain", 512, "LegacyCompact", False, 600, 3800),
        record("chain", 512, "OperatorPcg", True, 40, 600),
        record("chain", 512, "FixedNeighbor", True, 100, 400),
    ]
    records[0]["formal_result"]["status"] = "rss_limit"
    for item in records:
        item['input_sha256'] = 'same-input'
    analysis = frontier.analyze(records)
    group = analysis["groups"][0]
    assert group["fastest_route"] == "OperatorPcg"
    assert group["least_rss_route"] == "FixedNeighbor"
    assert group["formal_envelope_pass_routes"] == ["FixedNeighbor", "OperatorPcg"]
    assert group["matched_input"] is True
    assert analysis["route_position"] == "FixedNeighbor-memory-advantaged"


def test_search_completion_does_not_claim_full_endpoint_qualification():
    report = {
        "execution": {"status": "completed", "runs": [{
            "kind": "measurement", "status": "completed", "process_wall_seconds": 20,
            "peak_rss_bytes": 100 * 1024**2,
            "result": {"qualified": None, "details": {
                "search_completed": True, "search_seconds": 19, "profile_evaluations": 6,
            }},
        }]},
        "resources": {"elapsed_seconds": 19},
        "problem": {"atoms": 8, "voxels": 200, "parameters": 24},
        "result": {},
    }
    summary = frontier.result_summary(report, "search-only", 600)
    assert summary["completed"] is True
    assert summary["endpoint_certified"] is None


def test_least_rss_uses_formal_timeout_observation_when_diagnostic_is_deferred():
    records = [
        record("cube", 2048, "OperatorPcg", False, 600, 3200),
        record("cube", 2048, "FixedNeighbor", False, 600, 490),
    ]
    records[0]["diagnostic_result"] = {
        "status": "not-run", "completed": False,
        "reason": "extended diagnostic deferred until formal frontier is recorded",
    }
    records[1]["diagnostic_result"] = dict(records[1]["formal_result"])
    group = frontier.analyze(records)["groups"][0]
    assert group["fastest_route"] is None
    assert group["least_rss_route"] == "FixedNeighbor"


def test_2048_is_appended_only_after_all_lower_cases():
    assert set(frontier.DEFAULT_CASES) == {
        "chain-512", "cube-512", "chain-768", "cube-768", "chain-1024", "cube-1024",
    }
    for case in ("chain-2048", "cube-2048"):
        assert frontier._case(case)[1] == 2048


def test_output_writer_uses_both_resource_envelope_records():
    with tempfile.TemporaryDirectory() as temporary:
        manifest = {
            'formal_envelope': {'wall_seconds': 600},
            'diagnostic_envelope': {'wall_seconds': 7200},
        }
        frontier.write_outputs(Path(temporary), manifest, [])
        readme = (Path(temporary) / 'README.md').read_text()
        assert '600-second cap' in readme
        assert '7200-second wall cap' in readme


def test_skip_diagnostics_records_formal_timeouts_without_retrying():
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        build = root / 'build'
        build.mkdir()
        (build / 'CMakeCache.txt').write_text(
            'RHBM_GEM_JOINT_SPARSE_BACKEND:STRING=SPQR\n')
        output = root / 'frontier'
        args = frontier.build_parser().parse_args([
            '--build-dir', str(build), '--output-dir', str(output),
            '--cases', 'chain-512', '--search-only', '--skip-diagnostics',
        ])
        calls = []

        def formal_timeout(args, profile, case, route, timeout, report_path):
            calls.append((case, route, timeout))
            return {
                'rows': 512000, 'parameter_count': 1536,
                'result': {'status': 'timeout', 'completed': False,
                           'timeout_seconds': timeout, 'total_seconds': timeout},
            }

        with patch.object(frontier, '_run_case', formal_timeout):
            frontier.run_campaign(args)

        assert len(calls) == 3
        with (output / 'summary.csv').open(newline='') as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 3
        assert {row['diagnostic_status'] for row in rows} == {'not-run'}
        manifest = frontier.read(output / 'campaign-manifest.json')
        assert manifest['diagnostic_retry_policy'] == 'deferred'


if __name__ == '__main__':
    import sys
    import unittest

    tests = (
        test_frontier_reports_fastest_rss_and_formal_envelope_separately,
        test_search_completion_does_not_claim_full_endpoint_qualification,
        test_least_rss_uses_formal_timeout_observation_when_diagnostic_is_deferred,
        test_2048_is_appended_only_after_all_lower_cases,
        test_output_writer_uses_both_resource_envelope_records,
        test_skip_diagnostics_records_formal_timeouts_without_retrying,
    )
    suite = unittest.TestSuite(unittest.FunctionTestCase(test) for test in tests)
    sys.exit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
