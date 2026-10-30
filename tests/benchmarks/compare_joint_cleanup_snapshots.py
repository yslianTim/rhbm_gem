"""Compare cleanup snapshots without numerical tolerances (run in the same build environment)."""
import collections
import json
from pathlib import Path
import struct
import sys


def exact(left, right):
    if isinstance(left, float) or isinstance(right, float):
        return struct.pack("d", left) == struct.pack("d", right)
    if isinstance(left, list):
        return len(left) == len(right) and all(exact(a, b) for a, b in zip(left, right))
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(exact(left[k], right[k]) for k in left)
    return left == right


def compare(before, after):
    assert len(before) == len(after), "fixture count changed"
    totals = collections.Counter()
    for old, new in zip(before, after):
        for key in old:
            if key not in ("points", "conditioning", "wall_seconds", "cpu_seconds", "peak_rss_bytes"):
                assert exact(old[key], new[key]), (old["fixture"], key)
        assert len(old["conditioning"]) == len(new["conditioning"])
        for left, right in zip(old["conditioning"], new["conditioning"]):
            for key in left:
                if key == "workspace_peak_memory_estimate_bytes":
                    # Removed one double from EvaluationContext in each prepared block.
                    assert left[key] - right[key] == 8 * len(old["conditioning"])
                else:
                    assert exact(left[key], right[key]), (old["fixture"], key)
        transitions = collections.Counter()
        assert len(old["points"]) == len(new["points"])
        for left, right in zip(old["points"], new["points"]):
            for key in ("atom", "conditioning", "convergence", "validity"):
                assert exact(left[key], right[key]), (old["fixture"], key)
            if right["convergence"] != 0 or right["validity"] != 0:
                expected = "Unavailable"
            elif right["conditioning"] in ("Warning", "Unavailable"):
                expected = "Warning"
            elif right["conditioning"] == "NotAssessed":
                expected = "NotAssessed"
            else:
                assert right["conditioning"] == "Adequate"
                expected = "Available"
            assert right["quality"] == expected, (old["fixture"], right)
            transitions[left["quality"] + " -> " + right["quality"]] += 1
        totals.update(transitions)
        print(old["fixture"], "exact numerical state, conditioning and qualification;", dict(transitions))
    print("ALL", dict(totals))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: compare_joint_cleanup_snapshots.py <before.jsonl> <after.jsonl>")
    compare(*[[json.loads(line) for line in Path(path).read_text().splitlines()] for path in sys.argv[1:]])
