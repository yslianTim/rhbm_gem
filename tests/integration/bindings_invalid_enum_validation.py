import rhbm_gem_module as m


def has_parse_issue(command, option_name: str) -> bool:
    return any(
        issue.option_name == option_name and issue.phase == m.ValidationPhase.Parse
        for issue in command.GetValidationIssues()
    )


def main() -> int:

    dump = m.ResultDumpCommand()
    dump.SetPrinterChoice(999)
    dump.SetModelKeyTagList("model")
    assert not dump.PrepareForExecution()
    assert has_parse_issue(dump, "--printer")

    display = m.PotentialDisplayCommand()
    display.SetPainterChoice(999)
    display.SetModelKeyTagList("model")
    assert not display.PrepareForExecution()
    assert has_parse_issue(display, "--painter")

    model_test = m.HRLModelTestCommand()
    model_test.SetTesterChoice(999)
    assert not model_test.PrepareForExecution()
    assert has_parse_issue(model_test, "--tester")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
