#!/usr/bin/env python3
"""Generate a command scaffold.

By default this script only creates command/binding/test/doc skeleton files.
When ``--wire`` is set, it also updates the manifest and CMake source lists.
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
    command_id: str
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


def _append_command_entry(text: str, spec: ScaffoldSpec) -> tuple[str, bool]:
    if f"{spec.command_id}," in text:
        return text, False
    block = (
        "\n"
        "RHBM_GEM_COMMAND(\n"
        f"    {spec.command_id},\n"
        f'    "{spec.cli_name}",\n'
        f'    "{spec.description}",\n'
        f"    {spec.profile})\n"
    )
    return text.rstrip() + block, True


def _append_cmake_list_entry(text: str, variable_name: str, entry: str) -> tuple[str, bool]:
    pattern = re.compile(
        rf"(^set\({re.escape(variable_name)}\n)(.*?)(^\))",
        re.MULTILINE | re.DOTALL,
    )
    match = pattern.search(text)
    if not match:
        raise RuntimeError(f"could not find set({variable_name} ...) block")

    entry_line = f"    {entry}"
    body_lines = [line for line in match.group(2).splitlines() if line.strip()]
    if entry_line in body_lines:
        return text, False

    body_lines.append(entry_line)
    body_lines.sort(key=lambda line: line.strip().lower())
    updated_body = "\n".join(body_lines) + "\n"
    updated_text = text[:match.start(2)] + updated_body + text[match.start(3):]
    return updated_text, True


def _append_registry_include(text: str, spec: ScaffoldSpec) -> tuple[str, bool]:
    include_line = f'#include "{spec.command_type}.hpp"'
    if include_line in text:
        return text, False

    include_pattern = re.compile(r'^(#include\s+"[^"]+Command\.hpp")\n', re.MULTILINE)
    matches = list(include_pattern.finditer(text))
    if not matches:
        raise RuntimeError("could not find any command include in CommandRegistry.hpp")

    insert_at = matches[-1].end()
    updated = text[:insert_at] + include_line + "\n" + text[insert_at:]
    return updated, True


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

#include "CommandBase.hpp"

namespace rhbm_gem {{

struct {spec.command_type.removesuffix("Command")}Request;

class {spec.command_type}
    : public CommandWithRequest<{spec.command_type.removesuffix("Command")}Request>
{{
public:
    explicit {spec.command_type}(CommonOptionProfile profile);
    ~{spec.command_type}() override = default;

private:
    void NormalizeRequest() override;
    void ValidateOptions() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
}};

}} // namespace rhbm_gem
"""


def _source_template(spec: ScaffoldSpec) -> str:
    return f"""#include "internal/command/{spec.command_type}.hpp"

#include <rhbm_gem/core/command/CommandApi.hpp>

namespace rhbm_gem {{

{spec.command_type}::{spec.command_type}(CommonOptionProfile profile) :
    CommandWithRequest<{spec.command_type.removesuffix("Command")}Request>{{ profile }}
{{
}}

void {spec.command_type}::NormalizeRequest()
{{
    auto & request{{ MutableRequest() }};
    (void)request;
    // Normalize typed request fields here. Use CommandBase helpers when you need
    // validation issues, filesystem checks, or fallback behavior.
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

1. Add `{spec.command_type.removesuffix("Command")}` into `include/rhbm_gem/core/command/CommandList.def`.
   If it is experimental, place it inside the `RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` block.
2. Add the request struct in `include/rhbm_gem/core/command/CommandApi.hpp`.
   Define `VisitFields(...)` there so CLI and Python bindings pick it up automatically.
3. Add `#include "{spec.command_type}.hpp"` to `src/core/internal/command/CommandRegistry.hpp`.
"""


def _wire_registration_files(root: Path, spec: ScaffoldSpec, dry_run: bool, strict: bool) -> None:
    command_def = root / "include" / "rhbm_gem" / "core" / "command" / "CommandList.def"
    command_registry = root / "src" / "core" / "internal" / "command" / "CommandRegistry.hpp"
    source_cmake = root / "src" / "CMakeLists.txt"
    tests_cmake = root / "tests" / "CMakeLists.txt"

    _update_file(
        command_def,
        lambda text: _append_command_entry(text, spec),
        dry_run,
        strict,
        "Append a new RHBM_GEM_COMMAND(...) block to CommandList.def.",
    )
    _update_file(
        command_registry,
        lambda text: _append_registry_include(text, spec),
        dry_run,
        strict,
        f'Add #include "{spec.command_type}.hpp" to CommandRegistry.hpp.',
    )
    _update_file(
        source_cmake,
        lambda text: _append_cmake_list_entry(
            text,
            "RHBM_GEM_COMMAND_SOURCES",
            f"core/command/{spec.command_type}.cpp",
        ),
        dry_run,
        strict,
        f"Add core/command/{spec.command_type}.cpp to RHBM_GEM_COMMAND_SOURCES.",
    )
    _update_file(
        tests_cmake,
        lambda text: _append_cmake_list_entry(
            text,
            "RHBM_GEM_CORE_COMMAND_TEST_SOURCES",
            f"core/command/{spec.command_type}_test.cpp",
        ),
        dry_run,
        strict,
        f"Add core/command/{spec.command_type}_test.cpp to RHBM_GEM_CORE_COMMAND_TEST_SOURCES.",
    )


def build_spec(args: argparse.Namespace) -> ScaffoldSpec:
    base = _to_title_case(args.name)
    command_type = base if base.endswith("Command") else f"{base}Command"
    command_id = command_type.removesuffix("Command")
    cli_name = _to_cli_name(args.cli_name or base)
    description = args.description or f"Run {cli_name}"
    return ScaffoldSpec(
        command_id=command_id,
        command_type=command_type,
        cli_name=cli_name,
        description=description,
        profile=args.profile,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate command scaffold files.")
    parser.add_argument("--name", required=True, help="Base command name, e.g. Example or ExampleCommand.")
    parser.add_argument("--cli-name", help="CLI subcommand token. Defaults to name converted to snake_case.")
    parser.add_argument("--description", help="Command description text.")
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
        help=("Also update CommandList.def and the command/test CMake source lists."),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=("Only valid with --wire. Fail-fast when manifest or CMake wiring fails."),
    )
    args = parser.parse_args()
    if args.strict and not args.wire:
        parser.error("--strict requires --wire.")

    root = Path(__file__).resolve().parents[3]
    spec = build_spec(args)

    doc_stem = re.sub(r"(?<!^)([A-Z])", r"-\1", spec.command_type.removesuffix("Command")).lower()
    files = {
        root / "src" / "core" / "internal" / "command" / f"{spec.command_type}.hpp": _header_template(spec),
        root / "src" / "core" / "command" / f"{spec.command_type}.cpp": _source_template(spec),
        root / "tests" / "core" / "command" / f"{spec.command_type}_test.cpp": _test_template(spec),
        root / "docs" / "developer" / "commands" / f"{doc_stem}.md": _doc_template(spec),
    }

    for path, content in files.items():
        _write_new(path, content, args.dry_run)

    if args.wire:
        _wire_registration_files(root, spec, args.dry_run, args.strict)

    print("\nScaffold complete.")
    if args.wire:
        print("Registration/manifests and CMake source lists were wired automatically.")
    else:
        print("Next: wire CommandList.def, update CommandRegistry.hpp and CMake source lists, and add the public request schema.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
