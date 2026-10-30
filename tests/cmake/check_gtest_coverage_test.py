#!/usr/bin/env python3

"""Deterministic fixtures for the registration coverage guard; no test bodies run."""

from contextlib import redirect_stderr, redirect_stdout
import io
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import check_gtest_coverage as guard


ORDINARY = "Ordinary.\n  First\n  Second\n"
TYPED = "Typed/0.  # TypeParam = int\n  Maps\nTyped/1.  # TypeParam = float\n  Maps\n"
PARAMETERIZED = "Values/Param.\n  Maps/0  # GetParam() = 42\n"
FULL = "Running main() from fixture/gtest_main.cc\n" + ORDINARY + TYPED + PARAMETERIZED


class CoverageTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.binary = Path(self.directory.name) / "test binary"
        self.binary.touch()

    def audit(self, groups, full=FULL):
        entries = []
        outputs = {"--gtest_filter=*": full}
        for name, pattern, output in groups:
            argument = f"--gtest_filter={pattern}"
            entries.append({"name": name, "command": [str(self.binary), argument]})
            outputs[argument] = output

        def listing(command, **kwargs):
            self.assertIn("--gtest_list_tests", command)
            return outputs[next(arg for arg in command if arg.startswith("--gtest_filter="))]

        with patch.object(guard, "run", side_effect=listing), redirect_stdout(io.StringIO()):
            return guard.audit(self.binary, {"tests": entries})

    def report(self, result):
        output = io.StringIO()
        with redirect_stdout(output):
            status = guard.report(*result)
        return status, output.getvalue()

    def test_complete_including_typed_and_value_parameters(self):
        result = self.audit([("all", "Ordinary.*:Typed/*.*:*/Param.*", FULL)])
        self.assertEqual(len(result[0]), 5)
        self.assertEqual(self.report(result)[0], 0)

    def test_missing_ordinary(self):
        status, output = self.report(self.audit([("partial", "-Ordinary.*", TYPED + PARAMETERIZED)]))
        self.assertEqual(status, 1)
        self.assertIn("MISSING: Ordinary.First", output)
        self.assertIn("MISSING: Ordinary.Second", output)

    def test_missing_typed(self):
        status, output = self.report(self.audit([("old", "*Typed.*:Ordinary.*:*/Param.*", ORDINARY + PARAMETERIZED)]))
        self.assertEqual(status, 1)
        self.assertIn("MISSING: Typed/0.Maps", output)
        self.assertIn("MISSING: Typed/1.Maps", output)

    def test_duplicate_is_reported_without_failing(self):
        status, output = self.report(self.audit([
            ("all", "Ordinary.*:Typed/*.*:*/Param.*", FULL),
            ("overlap", "Ordinary.*", ORDINARY),
        ]))
        self.assertEqual(status, 0)
        self.assertIn("duplicate coverage: 2", output)
        self.assertIn("DUPLICATE (unexpected): Ordinary.First", output)

    def test_intentional_persistence_duplicate(self):
        output = self.report(({"DataObjectPersistenceTest.RoundTrip"}, set(), {
            "DataObjectPersistenceTest.RoundTrip": ["rhbm_tests_data_schema", "joint_persistence_contract_test"]
        }))[1]
        self.assertIn("DUPLICATE (intentional)", output)

    def test_disabled_are_explicitly_excluded(self):
        disabled = "DISABLED_Suite.\n  Case\nPrefix/DISABLED_Typed/0.\n  Case\nSuite.\n  DISABLED_Case/0\n"
        result = self.audit([("enabled", "Ordinary.*", ORDINARY)], full=ORDINARY + disabled)
        status, output = self.report(result)
        self.assertEqual(status, 0)
        self.assertEqual(len(result[1]), 3)
        self.assertIn("disabled excluded: 3", output)
        self.assertIn("not run by default", output)

    def test_bad_empty_or_incomplete_listings_fail(self):
        for listing in ("", "not a GoogleTest listing", "  Orphan", "Empty.\n", "Suite.\n  Case\nEmpty.\n", "Suite.\n  Case\n  Case\n"):
            with self.subTest(listing=listing), self.assertRaises(ValueError):
                guard.parse_tests(listing)

    def test_empty_filtered_listing_fails(self):
        with self.assertRaisesRegex(ValueError, "Empty or incomplete"):
            self.audit([("empty", "NoSuchSuite.*", "")])

    def test_executable_missing(self):
        self.binary.unlink()
        with self.assertRaisesRegex(ValueError, "does not exist"):
            guard.audit(self.binary, {"tests": []})

    def test_missing_or_invalid_filter(self):
        for arguments in ([], ["--gtest_filter="], ["--gtest_filter=*", "--gtest_filter=-*"]):
            with self.subTest(arguments=arguments), patch.object(guard, "run", return_value=FULL):
                with self.assertRaisesRegex(ValueError, "one explicit nonempty"):
                    guard.audit(self.binary, {"tests": [{"name": "bad", "command": [str(self.binary), *arguments]}]})

    def test_unavailable_registration_information(self):
        for entries in ({}, {"tests": []}, {"tests": [{"name": "missing command"}]},
                        {"tests": [{"name": "other", "command": ["/other/binary"]}]}):
            with self.subTest(entries=entries), patch.object(guard, "run", return_value=FULL):
                with self.assertRaises(ValueError):
                    guard.audit(self.binary, entries)

    def test_group_execution_failure_propagates(self):
        entry = {"name": "broken", "command": [str(self.binary), "--gtest_filter=Suite.*"]}
        with patch.object(guard, "run", side_effect=[FULL, ValueError("filter execution failed")]):
            with self.assertRaisesRegex(ValueError, "filter execution failed"):
                guard.audit(self.binary, {"tests": [entry]})

    def test_disabled_ctest_group_cannot_supply_coverage(self):
        entry = {"name": "disabled", "command": [str(self.binary), "--gtest_filter=*"],
                 "properties": [{"name": "DISABLED", "value": True}]}
        with patch.object(guard, "run", return_value=FULL), redirect_stdout(io.StringIO()):
            result = guard.audit(self.binary, {"tests": [entry]})
        self.assertEqual(self.report(result)[0], 1)
        self.assertEqual(result[2], {})

    def test_inconsistent_group_registry_fails(self):
        with self.assertRaisesRegex(ValueError, "outside the full registry"):
            self.audit([("wrong", "Other.*", "Other.\n  Case\n")])

    def test_invalid_ctest_json_fails_main(self):
        arguments = ["guard", "--binary", str(self.binary), "--build-dir", self.directory.name, "--ctest", "ctest"]
        with patch("sys.argv", arguments), patch.object(guard, "run", return_value="{incomplete"), redirect_stderr(io.StringIO()):
            self.assertEqual(guard.main(), 1)

    def test_process_errors_fail(self):
        for error in (OSError("unavailable"), subprocess.CalledProcessError(2, ["binary"], stderr="bad filter"),
                      subprocess.TimeoutExpired(["binary"], 30)):
            with self.subTest(error=error), patch.object(subprocess, "run", side_effect=error):
                with self.assertRaisesRegex(ValueError, "Cannot list registration"):
                    guard.run(["binary", "--gtest_list_tests"])

    def test_filter_and_environment_are_forwarded_without_shell(self):
        pattern = "--gtest_filter=Ordinary.*:Typed/*.*-Ordinary.Second"
        entry = {"name": "filtered", "command": [str(self.binary), pattern], "properties": [
            {"name": "ENVIRONMENT", "value": ["OMP_NUM_THREADS=1"]},
            {"name": "WORKING_DIRECTORY", "value": self.directory.name},
        ]}
        with patch.object(guard, "run", side_effect=[FULL, ORDINARY]), redirect_stdout(io.StringIO()):
            guard.audit(self.binary, {"tests": [entry]})
            command = guard.run.call_args.args[0]
            self.assertEqual(command, [str(self.binary), pattern, "--gtest_list_tests", "--gtest_color=no"])
            self.assertEqual(guard.run.call_args.kwargs["env"]["OMP_NUM_THREADS"], "1")
            self.assertEqual(guard.run.call_args.kwargs["cwd"], self.directory.name)


if __name__ == "__main__":
    unittest.main()
