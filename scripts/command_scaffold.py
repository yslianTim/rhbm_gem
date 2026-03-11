#!/usr/bin/env python3
"""Generate a built-in command scaffold.

By default this script only creates command/binding/test/doc skeleton files.
When ``--wire`` is set, it also updates registration/manifests in-place.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import sys
from typing import Callable


@dataclass(frozen=True)
class ScaffoldSpec:
    command_type: str
    cli_name: str
    description: str
    python_binding_name: str
    profile: str


@dataclass(frozen=True)
class WireTargets:
    builtins_def: Path
    catalog_cpp: Path
    source_manifest: Path
    bindings_cmake: Path
    binding_helpers_hpp: Path
    core_bindings_cpp: Path
    core_tests_cmake: Path


def _to_title_case(text: str) -> str:
    parts = re.split(r"[_\-\s]+", text.strip())
    return "".join(p[:1].upper() + p[1:] for p in parts if p)


def _to_cli_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")


def _write_new(path: Path, content: str, dry_run: bool) -> None:
    if path.exists():
        print(f"[skip] {path} already exists")
        return
    if dry_run:
        print(f"[new]  {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    print(f"[new]  {path}")


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _write_text(path: Path, content: str, dry_run: bool) -> None:
    if dry_run:
        print(f"[wire] {path}")
        return
    path.write_text(content, encoding="utf-8")
    print(f"[wire] {path}")


def _append_before_closing_paren_in_set(
    text: str,
    set_name: str,
    entry: str,
) -> tuple[str, bool]:
    if entry in text:
        return text, False
    lines = text.splitlines(keepends=True)
    in_target_set = False
    for i, line in enumerate(lines):
        if not in_target_set and re.match(rf"^\s*set\(\s*{re.escape(set_name)}\b", line):
            in_target_set = True
            continue
        if in_target_set and re.match(r"^\s*\)\s*$", line):
            lines.insert(i, f"    {entry}\n")
            return "".join(lines), True
    raise RuntimeError(f"Cannot find closing ')' for set({set_name} ...)")


def _insert_after_anchor_once(text: str, anchor: str, insertion: str) -> tuple[str, bool]:
    if insertion.strip() in text:
        return text, False
    idx = text.find(anchor)
    if idx < 0:
        raise RuntimeError(f"Anchor not found: {anchor}")
    insert_at = idx + len(anchor)
    return text[:insert_at] + insertion + text[insert_at:], True


def _append_builtins_entry(text: str, spec: ScaffoldSpec) -> tuple[str, bool]:
    if f"{spec.command_type}," in text:
        return text, False
    block = (
        "\n"
        "RHBM_GEM_BUILTIN_COMMAND(\n"
        f"    {spec.command_type},\n"
        f'    "{spec.cli_name}",\n'
        f'    "{spec.description}",\n'
        f'    "{spec.python_binding_name}")\n'
    )
    return text.rstrip() + block, True


def _ensure_line(text: str, line: str) -> tuple[str, bool]:
    if line in text:
        return text, False
    return text + ("\n" if not text.endswith("\n") else "") + line + "\n", True


def _update_file(
    path: Path,
    transform: Callable[[str], tuple[str, bool]],
    dry_run: bool,
) -> bool:
    original = _read_text(path)
    updated, changed = transform(original)
    if changed:
        _write_text(path, updated, dry_run)
    else:
        print(f"[skip] {path} already wired")
    return changed


def _header_template(spec: ScaffoldSpec) -> str:
    return f"""#pragma once

#include <rhbm_gem/core/command/CommandBase.hpp>

namespace CLI
{{
    class App;
}}

namespace rhbm_gem {{

struct {spec.command_type}Options : public CommandOptions
{{
}};

class {spec.command_type}
    : public CommandWithProfileOptions<
          {spec.command_type}Options,
          CommandId::{spec.command_type.removesuffix("Command")},
          CommonOptionProfile::{spec.profile}>
{{
public:
    explicit {spec.command_type}(const DataIoServices & data_io_services);
    ~{spec.command_type}() override = default;

private:
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    void ValidateOptions() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
}};

}} // namespace rhbm_gem
"""


def _source_template(spec: ScaffoldSpec) -> str:
    return f"""#include <rhbm_gem/core/command/{spec.command_type}.hpp>

namespace rhbm_gem {{

{spec.command_type}::{spec.command_type}(const DataIoServices & data_io_services) :
    CommandWithProfileOptions<
        {spec.command_type}Options,
        CommandId::{spec.command_type.removesuffix("Command")},
        CommonOptionProfile::{spec.profile}>{{ data_io_services }}
{{
}}

void {spec.command_type}::RegisterCLIOptionsExtend(CLI::App * cmd)
{{
    (void)cmd;
}}

void {spec.command_type}::ValidateOptions()
{{
}}

void {spec.command_type}::ResetRuntimeState()
{{
}}

bool {spec.command_type}::ExecuteImpl()
{{
    return true;
}}

}} // namespace rhbm_gem
"""


def _binding_template(spec: ScaffoldSpec) -> str:
    bind_fn = f'Bind{spec.command_type.removesuffix("Command")}'
    return f"""#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/{spec.command_type}.hpp>

namespace rhbm_gem::bindings {{

void {bind_fn}(py::module_ & module)
{{
    auto command{{ BindBuiltInCommand<{spec.command_type}>(module) }};
    command
        .def(py::init<const DataIoServices &>())
        .def("Execute", &{spec.command_type}::Execute);
    BindCommandDiagnostics(command);
}}

}} // namespace rhbm_gem::bindings
"""


def _test_template(spec: ScaffoldSpec) -> str:
    suite = f"{spec.command_type}Test"
    return f"""#include <gtest/gtest.h>

namespace {{

TEST({suite}, ScaffoldPlaceholder)
{{
    SUCCEED();
}}

}} // namespace
"""


def _doc_template(spec: ScaffoldSpec) -> str:
    return f"""# {spec.command_type}

Scaffold generated for CLI command `{spec.cli_name}`.

## Registration Checklist

1. Add `{spec.command_type}` into `src/core/internal/BuiltInCommandList.def`.
2. Include `{spec.command_type}.hpp` in `src/core/command/BuiltInCommandCatalog.cpp`.
3. Add source files to `src/rhbm_gem_sources.cmake` and `bindings/CMakeLists.txt`.
4. Wire binding function in `bindings/BindingHelpers.hpp` and `bindings/CoreBindings.cpp`.
5. Add tests to `tests/cmake/core_tests.cmake`.
"""


def _wire_registration_files(root: Path, spec: ScaffoldSpec, dry_run: bool) -> None:
    stem = spec.command_type.removesuffix("Command")
    bind_fn = f"Bind{stem}"
    binding_unit = f"{stem}Bindings.cpp"
    wire = WireTargets(
        builtins_def=root / "src" / "core" / "internal" / "BuiltInCommandList.def",
        catalog_cpp=root / "src" / "core" / "command" / "BuiltInCommandCatalog.cpp",
        source_manifest=root / "src" / "rhbm_gem_sources.cmake",
        bindings_cmake=root / "bindings" / "CMakeLists.txt",
        binding_helpers_hpp=root / "bindings" / "BindingHelpers.hpp",
        core_bindings_cpp=root / "bindings" / "CoreBindings.cpp",
        core_tests_cmake=root / "tests" / "cmake" / "core_tests.cmake",
    )

    _update_file(
        wire.builtins_def,
        lambda text: _append_builtins_entry(text, spec),
        dry_run,
    )
    _update_file(
        wire.catalog_cpp,
        lambda text: _insert_after_anchor_once(
            text,
            "#include <rhbm_gem/core/command/ResultDumpCommand.hpp>\n",
            f"#include <rhbm_gem/core/command/{spec.command_type}.hpp>\n",
        ),
        dry_run,
    )
    _update_file(
        wire.source_manifest,
        lambda text: _append_before_closing_paren_in_set(
            text,
            "RHBM_GEM_LIBRARY_SOURCES",
            f"core/command/{spec.command_type}.cpp",
        ),
        dry_run,
    )
    _update_file(
        wire.bindings_cmake,
        lambda text: _append_before_closing_paren_in_set(
            text,
            "BINDINGS_SOURCES",
            binding_unit,
        ),
        dry_run,
    )
    _update_file(
        wire.binding_helpers_hpp,
        lambda text: _insert_after_anchor_once(
            text,
            "void BindHRLModelTest(py::module_ & module);\n",
            f"void {bind_fn}(py::module_ & module);\n",
        ),
        dry_run,
    )
    _update_file(
        wire.core_bindings_cpp,
        lambda text: _insert_after_anchor_once(
            text,
            "    rhbm_gem::bindings::BindHRLModelTest(module);\n",
            f"    rhbm_gem::bindings::{bind_fn}(module);\n",
        ),
        dry_run,
    )
    _update_file(
        wire.core_tests_cmake,
        lambda text: _append_before_closing_paren_in_set(
            text,
            "CORE_COMMAND_TEST_SOURCES",
            f"core/command/{spec.command_type}_test.cpp",
        ),
        dry_run,
    )


def build_spec(args: argparse.Namespace) -> ScaffoldSpec:
    base = _to_title_case(args.name)
    command_type = base if base.endswith("Command") else f"{base}Command"
    cli_name = _to_cli_name(args.cli_name or base)
    description = args.description or f"Run {cli_name}"
    python_binding_name = args.python_binding_name or command_type
    return ScaffoldSpec(
        command_type=command_type,
        cli_name=cli_name,
        description=description,
        python_binding_name=python_binding_name,
        profile=args.profile,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate built-in command scaffold files.")
    parser.add_argument("--name", required=True, help="Base command name, e.g. Example or ExampleCommand.")
    parser.add_argument("--cli-name", help="CLI subcommand token. Defaults to name converted to snake_case.")
    parser.add_argument("--description", help="Built-in description text.")
    parser.add_argument("--python-binding-name", help="Python class name for BuiltInCommandList.def.")
    parser.add_argument(
        "--profile",
        default="FileWorkflow",
        choices=("FileWorkflow", "DatabaseWorkflow"),
        help="CommonOptionProfile for the generated command.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print planned files without writing.")
    parser.add_argument(
        "--wire",
        action="store_true",
        help=(
            "Also update command registration/manifests "
            "(BuiltInCommandList.def, CMake lists, bindings, and core test list)."
        ),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    spec = build_spec(args)

    bind_unit = f'{spec.command_type.removesuffix("Command")}Bindings.cpp'
    doc_stem = re.sub(r"(?<!^)([A-Z])", r"-\1", spec.command_type.removesuffix("Command")).lower()
    files = {
        root / "include" / "rhbm_gem" / "core" / "command" / f"{spec.command_type}.hpp": _header_template(spec),
        root / "src" / "core" / "command" / f"{spec.command_type}.cpp": _source_template(spec),
        root / "bindings" / bind_unit: _binding_template(spec),
        root / "tests" / "core" / "command" / f"{spec.command_type}_test.cpp": _test_template(spec),
        root / "docs" / "developer" / "commands" / f"{doc_stem}.md": _doc_template(spec),
    }

    for path, content in files.items():
        _write_new(path, content, args.dry_run)

    if args.wire:
        _wire_registration_files(root, spec, args.dry_run)

    print("\nScaffold complete.")
    if args.wire:
        print("Registration/manifests were wired automatically.")
    else:
        print("Next: wire BuiltInCommandList.def, CMake manifests, binding registration, and tests/cmake.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
