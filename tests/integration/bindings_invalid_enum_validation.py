import rhbm_gem_module as m


def has_parse_issue(report, option_name: str) -> bool:
    return any(
        issue.option_name == option_name and issue.phase == m.ValidationPhase.Parse
        for issue in report.validation_issues
    )


def main() -> int:
    dump = m.ResultDumpRequest()
    dump.model_key_tag_list = ""
    report = m.RunResultDump(dump)
    assert not report.prepared
    assert has_parse_issue(report, "--model-keylist")

    analysis = m.PotentialAnalysisRequest()
    analysis.saved_key_tag = ""
    report = m.RunPotentialAnalysis(analysis)
    assert not report.prepared
    assert has_parse_issue(report, "--save-key")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
