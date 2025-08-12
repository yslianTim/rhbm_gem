#pragma once

#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

class HRLModelTestCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
    {
        double fit_range_min{ 0.0 };
        double fit_range_max{ 1.0 };
        double alpha_r{ 0.1 };
        double alpha_g{ 0.2 };
    };

private:
    Options m_options{};

public:
    HRLModelTestCommand(void);
    ~HRLModelTestCommand() = default;
    bool Execute(void) override;
    bool ValidateOptions(void) const override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    const CommandOptions & GetOptions(void) const override { return m_options; }
    CommandOptions & GetOptions(void) override { return m_options; }

    void SetFitRangeMinimum(double value) { m_options.fit_range_min = value; }
    void SetFitRangeMaximum(double value) { m_options.fit_range_max = value; }
    void SetAlphaR(double value) { m_options.alpha_r = value; }
    void SetAlphaG(double value) { m_options.alpha_g = value; }
    void SetThreadSize(int value) { m_options.thread_size = value; }

private:

};
