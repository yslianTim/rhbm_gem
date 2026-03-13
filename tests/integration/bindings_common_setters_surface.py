import rhbm_gem_module as m

REQUEST_TYPES = (
    m.PotentialAnalysisRequest,
    m.PotentialDisplayRequest,
    m.ResultDumpRequest,
    m.MapSimulationRequest,
    m.MapVisualizationRequest,
    m.PositionEstimationRequest,
    m.HRLModelTestRequest,
)

EXPECTED_COMMON_FIELDS = {
    "thread_size",
    "verbose_level",
    "database_path",
    "folder_path",
}


def main() -> int:
    for request_type in REQUEST_TYPES:
        request = request_type()
        common = request.common
        missing = [field for field in EXPECTED_COMMON_FIELDS if not hasattr(common, field)]
        assert not missing, f"{request_type.__name__}.common missing fields: {missing}"
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
