#!/usr/bin/env python3
"""Quickstart for the rhbm_gem_module Python bindings."""

from __future__ import annotations

import inspect

import rhbm_gem_module as rgm


def print_enum_members(enum_name: str) -> None:
    enum_cls = getattr(rgm, enum_name)
    member_names = [name for name in dir(enum_cls) if name.isupper()]
    print(f"[enum] {enum_name}")
    for name in member_names:
        member = getattr(enum_cls, name)
        try:
            value = int(member)
        except (TypeError, ValueError):
            value = member
        print(f"  - {name} = {value}")


def print_available_setters(obj: object) -> None:
    method_names = [
        name for name, member in inspect.getmembers(obj, predicate=callable)
        if name.startswith("Set")
    ]
    print(f"[class] {obj.__class__.__name__}")
    if not method_names:
        print("  - (no Set* methods exposed)")
        return
    for name in method_names:
        print(f"  - {name}")


def main() -> int:
    print("Imported module: rhbm_gem_module")
    print("")
    print("Available enums:")
    for enum_name in (
        "PotentialModel",
        "PartialCharge",
        "PainterType",
        "PrinterType",
        "TesterType",
    ):
        print_enum_members(enum_name)
    print("")

    commands = [
        rgm.HRLModelTestCommand(),
        rgm.MapSimulationCommand(),
        rgm.MapVisualizationCommand(),
        rgm.PotentialAnalysisCommand(),
        rgm.PositionEstimationCommand(),
        rgm.ResultDumpCommand(),
        rgm.PotentialDisplayCommand(),
    ]

    print("Available command setters:")
    for command in commands:
        print_available_setters(command)
    print("")
    print("Next step:")
    print("  python3 examples/python/01_end_to_end_from_test_data.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
