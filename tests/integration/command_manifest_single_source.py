from pathlib import Path
import sys


def main() -> int:
    project_root = Path(__file__).resolve().parents[2]
    scripts_dir = project_root / "scripts"
    sys.path.insert(0, str(scripts_dir))

    import command_manifest
    import generate_command_artifacts
    import check_command_sync

    command_def = project_root / "src" / "core" / "internal" / "CommandList.def"
    entries_from_manifest = command_manifest.parse_commands(command_def)
    entries_from_generator = generate_command_artifacts.parse_commands(command_def)
    entries_from_sync = check_command_sync.parse_commands(command_def)

    assert entries_from_manifest, "Expected at least one command entry."
    assert entries_from_manifest == entries_from_generator
    assert entries_from_manifest == entries_from_sync
    assert generate_command_artifacts.parse_commands is command_manifest.parse_commands
    assert check_command_sync.parse_commands is command_manifest.parse_commands

    outputs = generate_command_artifacts.compute_expected_outputs(
        project_root,
        entries_from_manifest,
    )
    drift_paths = command_manifest.collect_drift_paths(outputs)
    assert isinstance(drift_paths, list)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
