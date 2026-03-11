import rhbm_gem_module as m


# Snapshot guards against accidental API removals.
EXPECTED_METHODS = {
    "PotentialAnalysisCommand": {
        "Execute",
        "PrepareForExecution",
        "HasValidationErrors",
        "GetValidationIssues",
        "SetAsymmetryFlag",
        "SetFitRangeMinimum",
        "SetFitRangeMaximum",
        "SetAlphaR",
        "SetAlphaG",
        "SetModelFilePath",
        "SetMapFilePath",
        "SetSavedKeyTag",
        "SetTrainingReportDir",
        "SetSamplingSize",
        "SetSamplingRangeMinimum",
        "SetSamplingRangeMaximum",
        "SetSimulationFlag",
        "SetSimulatedMapResolution",
    },
    "PotentialDisplayCommand": {
        "Execute",
        "PrepareForExecution",
        "HasValidationErrors",
        "GetValidationIssues",
        "SetPainterChoice",
        "SetModelKeyTagList",
        "SetRefModelKeyTagListMap",
        "SetPickChainID",
        "SetPickResidueType",
        "SetPickElementType",
        "SetVetoChainID",
        "SetVetoResidueType",
        "SetVetoElementType",
    },
    "ResultDumpCommand": {
        "Execute",
        "PrepareForExecution",
        "HasValidationErrors",
        "GetValidationIssues",
        "SetPrinterChoice",
        "SetModelKeyTagList",
        "SetMapFilePath",
    },
    "MapSimulationCommand": {
        "Execute",
        "PrepareForExecution",
        "HasValidationErrors",
        "GetValidationIssues",
        "SetModelFilePath",
        "SetPotentialModelChoice",
        "SetPartialChargeChoice",
        "SetCutoffDistance",
        "SetGridSpacing",
        "SetBlurringWidthList",
    },
    "MapVisualizationCommand": {
        "Execute",
        "PrepareForExecution",
        "HasValidationErrors",
        "GetValidationIssues",
        "SetModelFilePath",
        "SetMapFilePath",
        "SetAtomSerialID",
        "SetSamplingSize",
        "SetWindowSize",
    },
    "PositionEstimationCommand": {
        "Execute",
        "PrepareForExecution",
        "HasValidationErrors",
        "GetValidationIssues",
        "SetMapFilePath",
        "SetIterationCount",
        "SetKNNSize",
        "SetAlpha",
        "SetThresholdRatio",
        "SetDedupTolerance",
    },
    "HRLModelTestCommand": {
        "Execute",
        "PrepareForExecution",
        "HasValidationErrors",
        "GetValidationIssues",
        "SetTesterChoice",
        "SetFitRangeMinimum",
        "SetFitRangeMaximum",
        "SetAlphaR",
        "SetAlphaG",
    },
}


def main() -> int:
    services = m.DataIoServices()
    for class_name, expected in EXPECTED_METHODS.items():
        command = getattr(m, class_name)(services)
        missing = [method for method in expected if not hasattr(command, method)]
        assert not missing, f"{class_name} missing expected methods: {missing}"
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
