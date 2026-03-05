import rhbm_gem_module as m


def main() -> int:
    binding_names = list(m._built_in_command_names)

    assert len(binding_names) == 7
    assert hasattr(m, "LogLevel")
    assert hasattr(m, "ValidationPhase")
    assert hasattr(m, "ValidationIssue")
    assert hasattr(m, "PrinterType")
    assert hasattr(m.PrinterType, "ATOM_OUTLIER")
    assert hasattr(m, "TesterType")
    assert hasattr(m.TesterType, "BENCHMARK")
    assert all(hasattr(m, name) for name in binding_names)

    instances = [getattr(m, name)() for name in binding_names]
    assert all(
        all(
            hasattr(instance, method)
            for method in ("PrepareForExecution", "HasValidationErrors", "GetValidationIssues")
        )
        for instance in instances
    )

    map_visualization = m.MapVisualizationCommand()
    assert hasattr(map_visualization, "SetModelFilePath")
    assert hasattr(map_visualization, "SetMapFilePath")
    assert hasattr(map_visualization, "SetAtomSerialID")
    assert hasattr(map_visualization, "SetSamplingSize")
    assert hasattr(map_visualization, "SetWindowSize")

    position_estimation = m.PositionEstimationCommand()
    assert hasattr(position_estimation, "SetMapFilePath")
    assert hasattr(position_estimation, "SetIterationCount")
    assert hasattr(position_estimation, "SetKNNSize")
    assert hasattr(position_estimation, "SetAlpha")
    assert hasattr(position_estimation, "SetThresholdRatio")
    assert hasattr(position_estimation, "SetDedupTolerance")

    model_test = m.HRLModelTestCommand()
    assert hasattr(model_test, "SetTesterChoice")
    assert hasattr(model_test, "SetFitRangeMinimum")
    assert hasattr(model_test, "SetFitRangeMaximum")
    assert hasattr(model_test, "SetAlphaR")
    assert hasattr(model_test, "SetAlphaG")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
