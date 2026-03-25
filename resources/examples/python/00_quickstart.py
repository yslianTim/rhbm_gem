#!/usr/bin/env python3
"""Quickstart for the rhbm_gem_module Python bindings."""

from __future__ import annotations

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


def print_request_fields(request_name: str) -> None:
    request = getattr(rgm, request_name)()
    fields = [name for name in dir(request) if not name.startswith("_")]
    print(f"[request] {request_name}")
    for name in fields:
        if callable(getattr(request, name)):
            continue
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

    print("Available run functions:")
    for name in (
        "RunPotentialAnalysis",
        "RunPotentialDisplay",
        "RunResultDump",
        "RunMapSimulation",
        "RunMapVisualization",
        "RunPositionEstimation",
        "RunHRLModelTest",
    ):
        print(f"  - {name}")
    print("")

    print("Available request fields:")
    for request_name in (
        "MapSimulationRequest",
        "PotentialAnalysisRequest",
        "ResultDumpRequest",
    ):
        print_request_fields(request_name)
    print("")

    print("Next step:")
    print("  python3 resources/examples/python/01_end_to_end_from_test_data.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
