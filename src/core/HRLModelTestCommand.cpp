#include "HRLModelTestCommand.hpp"
#include "HRLModelHelper.hpp"
#include "ROOTHelper.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

#include <random>

namespace {
CommandRegistrar<HRLModelTestCommand> registrar_model_test{
    "model_test",
    "Run HRL model simulation test"};
}

HRLModelTestCommand::HRLModelTestCommand(void) :
    CommandBase()
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::HRLModelTestCommand() called.");
}

void HRLModelTestCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::RegisterCLIOptionsExtend() called.");
    cmd->add_option("--fit-min", m_options.fit_range_min,
        "Minimum fitting range")->default_val(m_options.fit_range_min);
    cmd->add_option("--fit-max", m_options.fit_range_max,
        "Maximum fitting range")->default_val(m_options.fit_range_max);
    cmd->add_option("--alpha-r", m_options.alpha_r,
        "Alpha value for R")->default_val(m_options.alpha_r);
    cmd->add_option("--alpha-g", m_options.alpha_g,
        "Alpha value for G")->default_val(m_options.alpha_g);
}

bool HRLModelTestCommand::Execute(void)
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::Execute() called.");

    return true;
}

bool HRLModelTestCommand::ValidateOptions(void) const
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::ValidateOptions() called");
    if (m_options.fit_range_min >= m_options.fit_range_max)
    {
        Logger::Log(LogLevel::Error, "Invalid fitting range");
        return false;
    }
    return true;
}
