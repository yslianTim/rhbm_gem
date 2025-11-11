#pragma once

#include <CLI/CLI.hpp>
#include <filesystem>

#include "CommandBase.hpp"
#include "OptionEnumClass.hpp"

class HRLModelTestCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
    {
        TesterType tester_choice{ TesterType::BENCHMARK };
        double fit_range_min{ 0.0 };
        double fit_range_max{ 1.0 };
        double alpha_r{ 0.1 };
        double alpha_g{ 0.2 };
    };

private:
    Options m_options;

public:
    HRLModelTestCommand(void);
    ~HRLModelTestCommand() = default;
    bool Execute(void) override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    const CommandOptions & GetOptions(void) const override { return m_options; }
    CommandOptions & GetOptions(void) override { return m_options; }

    void SetTesterChoice(TesterType value);
    void SetFitRangeMinimum(double value);
    void SetFitRangeMaximum(double value);
    void SetAlphaR(double value);
    void SetAlphaG(double value);

private:
    void RunSimulationTestOnBenchMark(void);
    void RunSimulationTestOnDataOutlier(void);
    void RunSimulationTestOnMemberOutlier(void);
    void RunSimulationTestOnModelAlphaData(void);
    void RunSimulationTestOnModelAlphaMember(void);

};
