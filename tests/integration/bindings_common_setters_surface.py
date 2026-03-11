import rhbm_gem_module as m


# Derived from command profiles: all built-ins expose threading/verbose/output-folder,
# and database commands additionally expose database-path setter.
EXPECTED_COMMON_SETTERS = {
    "PotentialAnalysisCommand": {"SetThreadSize", "SetVerboseLevel", "SetFolderPath", "SetDatabasePath"},
    "PotentialDisplayCommand": {"SetThreadSize", "SetVerboseLevel", "SetFolderPath", "SetDatabasePath"},
    "ResultDumpCommand": {"SetThreadSize", "SetVerboseLevel", "SetFolderPath", "SetDatabasePath"},
    "MapSimulationCommand": {"SetThreadSize", "SetVerboseLevel", "SetFolderPath"},
    "MapVisualizationCommand": {"SetThreadSize", "SetVerboseLevel", "SetFolderPath"},
    "PositionEstimationCommand": {"SetThreadSize", "SetVerboseLevel", "SetFolderPath"},
    "HRLModelTestCommand": {"SetThreadSize", "SetVerboseLevel", "SetFolderPath"},
}


def main() -> int:
    services = m.DataIoServices()
    for class_name, expected_setters in EXPECTED_COMMON_SETTERS.items():
        command = getattr(m, class_name)(services)
        missing = [setter for setter in expected_setters if not hasattr(command, setter)]
        assert not missing, f"{class_name} missing common setters: {missing}"
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
