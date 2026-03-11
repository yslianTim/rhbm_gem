#!/usr/bin/env python3
"""Generate a built-in command scaffold.

By default this script only creates command/binding/test/doc skeleton files.
When ``--wire`` is set, it also updates registration/manifests in-place.
Use ``--wire --strict`` to fail-fast on anchor drift with repair hints.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import subprocess
import sys
from typing import Callable


@dataclass(frozen=True)
class ScaffoldSpec:
    command_type: str
    cli_name: str
    description: str
    profile: str


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


def _append_builtins_entry(text: str, spec: ScaffoldSpec) -> tuple[str, bool]:
    if f"{spec.command_type}," in text:
        return text, False
    block = (
        "\n"
        "RHBM_GEM_BUILTIN_COMMAND(\n"
        f"    {spec.command_type},\n"
        f'    "{spec.cli_name}",\n'
        f'    "{spec.description}")\n'
    )
    return text.rstrip() + block, True


def _update_file(
    path: Path,
    transform: Callable[[str], tuple[str, bool]],
    dry_run: bool,
    strict: bool = False,
    suggestion: str = "",
) -> bool:
    original = _read_text(path)
    try:
        updated, changed = transform(original)
    except RuntimeError as exc:
        hint = f" Suggested fix: {suggestion}" if suggestion else ""
        message = f"{path}: {exc}.{hint}"
        if strict:
            raise RuntimeError(message) from exc
        print(f"[warn] {message}")
        return False
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
    return f"""#include "BindingHelpers.hpp"

#include <rhbm_gem/core/command/{spec.command_type}.hpp>

namespace rhbm_gem::bindings {{

template <>
void BindCommand<{spec.command_type}>(py::module_ & module)
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
2. Run `python3 scripts/generate_builtin_command_artifacts.py`.
3. Keep generated artifacts clean (`python3 scripts/check_builtin_command_sync.py`).
"""


def _wire_registration_files(root: Path, spec: ScaffoldSpec, dry_run: bool, strict: bool) -> None:
    builtins_def = root / "src" / "core" / "internal" / "BuiltInCommandList.def"
    _update_file(
        builtins_def,
        lambda text: _append_builtins_entry(text, spec),
        dry_run,
        strict,
        "Append a new RHBM_GEM_BUILTIN_COMMAND(...) block to BuiltInCommandList.def.",
    )
    if dry_run:
        print("[wire] scripts/generate_builtin_command_artifacts.py")
        return

    generator = root / "scripts" / "generate_builtin_command_artifacts.py"
    result = subprocess.run(
        [sys.executable, str(generator)],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.stdout.strip():
        print(result.stdout.strip())
    if result.returncode != 0:
        if result.stderr.strip():
            print(result.stderr.strip(), file=sys.stderr)
        if strict:
            raise RuntimeError("Generator failed while wiring built-in command artifacts.")


def build_spec(args: argparse.Namespace) -> ScaffoldSpec:
    base = _to_title_case(args.name)
    command_type = base if base.endswith("Command") else f"{base}Command"
    cli_name = _to_cli_name(args.cli_name or base)
    description = args.description or f"Run {cli_name}"
    return ScaffoldSpec(
        command_type=command_type,
        cli_name=cli_name,
        description=description,
        profile=args.profile,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate built-in command scaffold files.")
    parser.add_argument("--name", required=True, help="Base command name, e.g. Example or ExampleCommand.")
    parser.add_argument("--cli-name", help="CLI subcommand token. Defaults to name converted to snake_case.")
    parser.add_argument("--description", help="Built-in description text.")
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
            "Also update BuiltInCommandList.def and regenerate derived built-in artifacts."
        ),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Only valid with --wire. Fail-fast when manifest generation fails."
        ),
    )
    args = parser.parse_args()
    if args.strict and not args.wire:
        parser.error("--strict requires --wire.")

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
        _wire_registration_files(root, spec, args.dry_run, args.strict)

    print("\nScaffold complete.")
    if args.wire:
        print("Registration/manifests were wired automatically (manifest + generated artifacts).")
    else:
        print("Next: wire BuiltInCommandList.def and run artifact generator.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
