#!/usr/bin/env python3
"""Validate built-in command manifest synchronization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import sys


@dataclass(frozen=True)
class CommandEntry:
    command_type: str
    cli_name: str
    description: str
    python_binding_name: str


def command_doc_stem(command_type: str) -> str:
    stem = command_type.removesuffix("Command")
    return re.sub(r"(?<!^)([A-Z])", r"-\1", stem).lower()


def parse_builtins(path: Path) -> list[CommandEntry]:
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(
        r"RHBM_GEM_BUILTIN_COMMAND\(\s*"
        r"(?P<command>[A-Za-z_][A-Za-z0-9_]*)\s*,\s*"
        r"\"(?P<cli>[^\"]+)\"\s*,\s*"
        r"\"(?P<desc>[^\"]*)\"\s*,\s*"
        r"\"(?P<py>[^\"]+)\"\s*"
        r"\)",
        re.MULTILINE,
    )
    return [
        CommandEntry(
            command_type=m.group("command"),
            cli_name=m.group("cli"),
            description=m.group("desc"),
            python_binding_name=m.group("py"),
        )
        for m in pattern.finditer(text)
    ]


def ensure_relative_order(text: str, tokens: list[str], label: str, errors: list[str]) -> None:
    positions: list[tuple[int, str]] = []
    for token in tokens:
        idx = text.find(token)
        if idx < 0:
            errors.append(f"{label} missing ordered token: {token}")
            return
        positions.append((idx, token))

    last = -1
    for idx, token in positions:
        if idx < last:
            ordered = " -> ".join(t for _, t in positions)
            errors.append(f"{label} has order drift (expected built-in order): {ordered}")
            return
        last = idx


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    builtins = parse_builtins(root / "src" / "core" / "internal" / "BuiltInCommandList.def")
    if not builtins:
        print("No built-in commands parsed from BuiltInCommandList.def")
        return 1

    catalog_text = (root / "src" / "core" / "command" / "BuiltInCommandCatalog.cpp").read_text(encoding="utf-8")
    sources_text = (root / "src" / "rhbm_gem_sources.cmake").read_text(encoding="utf-8")
    core_tests_text = (root / "tests" / "cmake" / "core_tests.cmake").read_text(encoding="utf-8")
    bindings_cmake_text = (root / "bindings" / "CMakeLists.txt").read_text(encoding="utf-8")
    binding_helpers_text = (root / "bindings" / "BindingHelpers.hpp").read_text(encoding="utf-8")
    core_bindings_text = (root / "bindings" / "CoreBindings.cpp").read_text(encoding="utf-8")

    errors: list[str] = []
    for entry in builtins:
        command = entry.command_type
        stem = command.removesuffix("Command")
        bind_fn = f"Bind{stem}"
        binding_unit = f"{stem}Bindings.cpp"
        expected_paths = [
            root / "include" / "rhbm_gem" / "core" / "command" / f"{command}.hpp",
            root / "src" / "core" / "command" / f"{command}.cpp",
            root / "bindings" / binding_unit,
            root / "tests" / "core" / "command" / f"{command}_test.cpp",
            root / "docs" / "developer" / "commands" / f"{command_doc_stem(command)}.md",
        ]
        for p in expected_paths:
            if not p.exists():
                errors.append(f"missing file for {command}: {p.relative_to(root)}")

        header_include = f"#include <rhbm_gem/core/command/{command}.hpp>"
        if header_include not in catalog_text:
            errors.append(f"BuiltInCommandCatalog.cpp missing include for {command}")

        source_ref = f"core/command/{command}.cpp"
        if source_ref not in sources_text:
            errors.append(f"rhbm_gem_sources.cmake missing source entry: {source_ref}")

        test_ref = f"core/command/{command}_test.cpp"
        if test_ref not in core_tests_text:
            errors.append(f"core_tests.cmake missing test entry: {test_ref}")

        if binding_unit not in bindings_cmake_text:
            errors.append(f"bindings/CMakeLists.txt missing source entry: {binding_unit}")

        helper_decl = f"void {bind_fn}(py::module_ & module);"
        if helper_decl not in binding_helpers_text:
            errors.append(f"BindingHelpers.hpp missing declaration: {helper_decl}")

        bind_call = f"rhbm_gem::bindings::{bind_fn}(module);"
        if bind_call not in core_bindings_text:
            errors.append(f"CoreBindings.cpp missing module registration call: {bind_call}")

        binding_text = (root / "bindings" / binding_unit).read_text(encoding="utf-8")
        bind_template_ref = f"BindBuiltInCommand<{command}>"
        if bind_template_ref not in binding_text:
            errors.append(f"{binding_unit} missing built-in bind template reference: {bind_template_ref}")

        doc_path = root / "docs" / "developer" / "commands" / f"{command_doc_stem(command)}.md"
        if doc_path.exists():
            doc_text = doc_path.read_text(encoding="utf-8")
            if f"# {command}" not in doc_text:
                errors.append(f"{doc_path.relative_to(root)} missing title header '# {command}'")
            if "## Registration Checklist" not in doc_text:
                errors.append(f"{doc_path.relative_to(root)} missing 'Registration Checklist' section")

    builtins_order = [entry.command_type for entry in builtins]
    builtins_stems = [command.removesuffix("Command") for command in builtins_order]
    ensure_relative_order(
        catalog_text,
        [f"#include <rhbm_gem/core/command/{command}.hpp>" for command in builtins_order],
        "BuiltInCommandCatalog.cpp include list",
        errors,
    )
    ensure_relative_order(
        binding_helpers_text,
        [f"void Bind{stem}(py::module_ & module);" for stem in builtins_stems],
        "BindingHelpers.hpp declarations",
        errors,
    )
    ensure_relative_order(
        core_bindings_text,
        [f"rhbm_gem::bindings::Bind{stem}(module);" for stem in builtins_stems],
        "CoreBindings.cpp bind calls",
        errors,
    )
    ensure_relative_order(
        bindings_cmake_text,
        [f"{stem}Bindings.cpp" for stem in builtins_stems],
        "bindings/CMakeLists.txt sources",
        errors,
    )
    ensure_relative_order(
        core_tests_text,
        [f"core/command/{command}_test.cpp" for command in builtins_order],
        "tests/cmake/core_tests.cmake command tests",
        errors,
    )

    if errors:
        print("Built-in command sync check failed:")
        for err in errors:
            print(f" - {err}")
        return 1

    print(f"Built-in command sync check passed ({len(builtins)} commands).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
