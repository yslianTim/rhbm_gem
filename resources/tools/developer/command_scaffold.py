#!/usr/bin/env python3
"""Generate a command scaffold.

By default this script creates command/binding/doc skeleton files.
When ``--wire`` is set, it also updates the command catalog and source CMake list.
Command tests live in grouped files under ``tests/core/command/`` and are updated manually.
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


def _is_experimental_command(spec: ScaffoldSpec) -> bool:
    return spec.command_type in {"MapVisualizationCommand", "PositionEstimationCommand"}


def _append_command_registry_entry(text: str, spec: ScaffoldSpec) -> tuple[str, bool]:
    if f"CommandEntry<{spec.command_id}Request>" in text:
        return text, False

    entry = (
        f"    CommandEntry<{spec.command_id}Request>{{\n"
        f'        "{spec.cli_name}",\n'
        f'        "{spec.description}",\n'
        f'        "{spec.command_id}Request",\n'
        f"        Execute{spec.command_type}\n"
        f"    }},\n"
    )

    if _is_experimental_command(spec):
        list_pattern = re.compile(
            r"(inline constexpr auto kExperimentalCommands = std::tuple\{\n)(.*?)(\};)",
            re.DOTALL,
        )
    else:
        list_pattern = re.compile(
            r"(inline constexpr auto kStableCommands = std::tuple\{\n)(.*?)(\};)",
            re.DOTALL,
        )

    list_match = list_pattern.search(text)
    if not list_match:
        command_group = "experimental" if _is_experimental_command(spec) else "stable"
        raise RuntimeError(f"could not find {command_group} command catalog in CommandCatalog.hpp")
    updated = text[:list_match.start(3)] + entry + text[list_match.start(3):]
    return updated, True


def _append_command_executor_declaration(text: str, spec: ScaffoldSpec) -> tuple[str, bool]:
    declaration = (
        f"CommandResult Execute{spec.command_type}"
        f"(const {spec.command_id}Request & request);"
    )
    if declaration in text:
        return text, False

    if _is_experimental_command(spec):
        pattern = re.compile(
            r"(#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE\n"
            r"CommandResult ExecuteMapVisualizationCommand.*?\n"
            r"CommandResult ExecutePositionEstimationCommand.*?\n)"
            r"(#endif)",
            re.DOTALL,
        )
    else:
        pattern = re.compile(
            r"(CommandResult ExecutePotentialAnalysisCommand.*?\n"
            r"CommandResult ExecutePotentialDisplayCommand.*?\n"
            r"CommandResult ExecuteResultDumpCommand.*?\n"
            r"CommandResult ExecuteMapSimulationCommand.*?\n"
            r"CommandResult ExecuteRHBMTestCommand.*?\n)"
            r"(#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE)",
            re.DOTALL,
        )

    match = pattern.search(text)
    if not match:
        command_group = "experimental" if _is_experimental_command(spec) else "stable"
        raise RuntimeError(f"could not find {command_group} executor declarations")
    updated = text[:match.start(2)] + declaration + "\n" + text[match.start(2):]
    return updated, True


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


def _source_template(spec: ScaffoldSpec) -> str:
    return f"""#include "detail/CommandBase.hpp"

namespace rhbm_gem::core {{

class {spec.command_type} final : public CommandBase<{spec.command_id}Request>
{{
public:
    {spec.command_type}();

private:
    void NormalizeAndValidateRequest({spec.command_id}Request & request) override;
    void ValidatePreparedRequest(const {spec.command_id}Request & request) override;
    bool ExecuteImpl(const {spec.command_id}Request & request) override;
}};

{spec.command_type}::{spec.command_type}() :
    CommandBase<{spec.command_id}Request>{{}}
{{
}}

void {spec.command_type}::NormalizeAndValidateRequest({spec.command_id}Request & request)
{{
    (void)request;
    // Normalize typed request fields and emit validation issues here.
}}

void {spec.command_type}::ValidatePreparedRequest(const {spec.command_id}Request & request)
{{
    (void)request;
}}

bool {spec.command_type}::ExecuteImpl(const {spec.command_id}Request & request)
{{
    (void)request;
    return true;
}}

namespace command_internal {{

CommandResult Execute{spec.command_type}(const {spec.command_id}Request & request)
{{
    {spec.command_type} command;
    return command.ExecuteRequest(request);
}}

}} // namespace command_internal

}} // namespace rhbm_gem::core
"""


def _doc_template(spec: ScaffoldSpec) -> str:
    return f"""# {spec.command_type}

Scaffold generated for CLI command `{spec.cli_name}`.

## Registration Checklist

1. Add the `{spec.command_id}Request` struct in `include/rhbm_gem/core/CommandTypes.hpp`.
   If it is experimental, keep it inside the `RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` block.
2. Add the internal request field catalog specialization in `src/core/command/detail/CommandCatalog.hpp`.
3. Add `Execute{spec.command_type}` declaration and `CommandEntry<{spec.command_id}Request>` to `src/core/command/detail/CommandCatalog.hpp`.
4. Add `core/command/{spec.command_type}.cpp` to the command source list in `src/CMakeLists.txt`.
5. Extend the grouped command tests under `tests/core/command/`.
   Command-specific validation, workflow, and output coverage usually belongs in
   `CommandScenarios_test.cpp`.
"""


def _wire_registration_files(root: Path, spec: ScaffoldSpec, dry_run: bool, strict: bool) -> None:
    command_catalog = root / "src" / "core" / "command" / "detail" / "CommandCatalog.hpp"
    source_cmake = root / "src" / "CMakeLists.txt"

    _update_file(
        command_catalog,
        lambda text: _append_command_executor_declaration(text, spec),
        dry_run,
        strict,
        f"Declare Execute{spec.command_type}(const {spec.command_id}Request &) in CommandCatalog.hpp.",
    )
    _update_file(
        command_catalog,
        lambda text: _append_command_registry_entry(text, spec),
        dry_run,
        strict,
        "Append a new typed CommandEntry block to CommandCatalog.hpp.",
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
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate command scaffold files.")
    parser.add_argument("--name", required=True, help="Base command name, e.g. Example or ExampleCommand.")
    parser.add_argument("--cli-name", help="CLI subcommand token. Defaults to name converted to snake_case.")
    parser.add_argument("--description", help="Command description text.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned files without writing.")
    parser.add_argument(
        "--wire",
        action="store_true",
        help=("Also update CommandCatalog.hpp and the command source CMake list."),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=("Only valid with --wire. Fail-fast when command-catalog or CMake wiring fails."),
    )
    args = parser.parse_args()
    if args.strict and not args.wire:
        parser.error("--strict requires --wire.")

    root = Path(__file__).resolve().parents[3]
    spec = build_spec(args)

    doc_stem = re.sub(r"(?<!^)([A-Z])", r"-\1", spec.command_type.removesuffix("Command")).lower()
    files = {
        root / "src" / "core" / "command" / f"{spec.command_type}.cpp": _source_template(spec),
        root / "docs" / "developer" / "commands" / f"{doc_stem}.md": _doc_template(spec),
    }

    for path, content in files.items():
        _write_new(path, content, args.dry_run)

    grouped_test_dir = root / "tests" / "core" / "command"
    print("\nTest guidance:")
    print(f"  - Extend grouped command tests under {grouped_test_dir}")
    print(
        "  - Command-specific validation, workflow, and output coverage usually belongs in "
        f"{grouped_test_dir / 'CommandScenarios_test.cpp'}"
    )

    if args.wire:
        _wire_registration_files(root, spec, args.dry_run, args.strict)

    print("\nScaffold complete.")
    if args.wire:
        print("Command catalog and command source CMake list were wired automatically.")
    else:
        print(
            "Next: wire CommandCatalog.hpp and source CMake list, "
            "add the public request DTO in CommandTypes.hpp and the CommandCatalog request fields, "
            "and extend the grouped command tests."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
