#!/usr/bin/env python3
"""Verify current convergence-audit authority and archive links."""

from __future__ import annotations

from pathlib import Path
import re
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEVELOPER_DOCS = PROJECT_ROOT / "docs" / "developer"
CURRENT_AUDIT = DEVELOPER_DOCS / "second-stage-outer-iteration-algorithm-audit.md"
HISTORY = DEVELOPER_DOCS / "audit-history" / "second-stage-convergence"
HISTORICAL_NAMES = (
    "convergence-safeguard-audit.md",
    "stationarity-active-coordinate-audit.md",
    "counterfactual-convergence-continuation-audit.md",
    "convergence-exposure-counterfactual-outcome-audit.md",
)
LINK_PATTERN = re.compile(r"]\(([^)]+)\)")


class ConvergenceAuditDocumentationTest(unittest.TestCase):
    def test_current_entry_points_use_one_authority(self) -> None:
        for path in (
            DEVELOPER_DOCS / "README.md",
            DEVELOPER_DOCS / "second-stage-local-fitting.md",
        ):
            text = path.read_text(encoding="utf-8")
            self.assertIn(CURRENT_AUDIT.name, text)
            for name in HISTORICAL_NAMES:
                self.assertNotIn(name, text)

    def test_historical_records_are_archived_with_front_matter(self) -> None:
        for name in HISTORICAL_NAMES:
            text = (HISTORY / name).read_text(encoding="utf-8")
            self.assertTrue(text.startswith(
                "---\n"
                "Status: Historical audit record.\n"
                "Current policy: second-stage-outer-iteration-algorithm-audit.md\n"
                "---\n"))
            self.assertIn(f"({HISTORY.relative_to(DEVELOPER_DOCS)}/{name})",
                          CURRENT_AUDIT.read_text(encoding="utf-8"))

    def test_local_markdown_links_resolve_after_move(self) -> None:
        paths = [CURRENT_AUDIT, *[HISTORY / name for name in HISTORICAL_NAMES]]
        for source in paths:
            for target in LINK_PATTERN.findall(source.read_text(encoding="utf-8")):
                target = target.split("#", 1)[0]
                if not target or "://" in target:
                    continue
                self.assertTrue(
                    (source.parent / target).resolve().is_file(),
                    f"{source}: unresolved link {target}")


if __name__ == "__main__":
    unittest.main()
