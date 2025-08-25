#include "HRLModelTestCommand.hpp"
#include "HRLModelHelper.hpp"
#include "ROOTHelper.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

#include <random>
#include <memory>
#include <vector>

namespace {
CommandRegistrar<HRLModelTestCommand> registrar_model_test{
    "model_test",
    "Run HRL model simulation test"};
}

HRLModelTestCommand::HRLModelTestCommand(void) :
    CommandBase(), m_options{}
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::HRLModelTestCommand() called.");
}

void HRLModelTestCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::RegisterCLIOptionsExtend() called.");
    std::map<std::string, TesterType> tester_map
    {
        {"1", TesterType::DATA_OUTLIER},      {"data_outlier",   TesterType::DATA_OUTLIER},
        {"2", TesterType::MEMBER_OUTLIER},    {"member_outlier", TesterType::MEMBER_OUTLIER},
        {"3", TesterType::MODEL_ALPHA_DATA},  {"alpha_data",     TesterType::MODEL_ALPHA_DATA},
        {"4", TesterType::MODEL_ALPHA_MEMBER},{"alpha_member",   TesterType::MODEL_ALPHA_MEMBER}
    };
    cmd->add_option_function<TesterType>("-t,--tester",
        [&](TesterType value) { SetTesterChoice(value); },
        "Tester option")
        ->default_val(TesterType::DATA_OUTLIER)
        ->transform(CLI::CheckedTransformer(tester_map, CLI::ignore_case));
    cmd->add_option_function<double>("--fit-min",
        [&](double value) { SetFitRangeMinimum(value); },
        "Minimum fitting range")->default_val(m_options.fit_range_min);
    cmd->add_option_function<double>("--fit-max",
        [&](double value) { SetFitRangeMaximum(value); },
        "Maximum fitting range")->default_val(m_options.fit_range_max);
    cmd->add_option_function<double>("--alpha-r",
        [&](double value) { SetAlphaR(value); },
        "Alpha value for R")->default_val(m_options.alpha_r);
    cmd->add_option_function<double>("--alpha-g",
        [&](double value) { SetAlphaG(value); },
        "Alpha value for G")->default_val(m_options.alpha_g);
}

bool HRLModelTestCommand::Execute(void)
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::Execute() called.");
    switch (m_options.tester_choice)
    {
        case TesterType::DATA_OUTLIER:
            RunSimulationTestOnDataOutlier();
            break;
        case TesterType::MEMBER_OUTLIER:
            RunSimulationTestOnMemberOutlier();
            break;
        case TesterType::MODEL_ALPHA_DATA:
            RunSimulationTestOnModelAlphaData();
            break;
        case TesterType::MODEL_ALPHA_MEMBER:
            RunSimulationTestOnModelAlphaMember();
            break;
        default:
            Logger::Log(LogLevel::Warning,
                        "Invalid tester choice input : ["
                        + std::to_string(static_cast<int>(m_options.tester_choice)) + "]");
            Logger::Log(LogLevel::Warning,
                        "Available Tester Choices:\n"
                        "  [1] Simulation Test on Data Outlier\n"
                        "  [2] Simulation Test on Member Outlier\n"
                        "  [3] Simulation Test on Model alpha_data\n"
                        "  [4] Simulation Test on Model alpha_member");
            break;
    }
    return true;
}

void HRLModelTestCommand::SetTesterChoice(TesterType value)
{
    m_options.tester_choice = value;
}

void HRLModelTestCommand::SetFitRangeMinimum(double value)
{
    m_options.fit_range_min = value;
}

void HRLModelTestCommand::SetFitRangeMaximum(double value)
{
    m_options.fit_range_max = value;
}

void HRLModelTestCommand::SetAlphaR(double value)
{
    m_options.alpha_r = value;
}

void HRLModelTestCommand::SetAlphaG(double value)
{
    m_options.alpha_g = value;
}

void HRLModelTestCommand::RunSimulationTestOnDataOutlier(void)
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::RunSimulationTestOnDataOutlier() called");
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnDataOutlier");

}

void HRLModelTestCommand::RunSimulationTestOnMemberOutlier(void)
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::RunSimulationTestOnMemberOutlier() called");
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnMemberOutlier");

}

void HRLModelTestCommand::RunSimulationTestOnModelAlphaData(void)
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::RunSimulationTestOnModelAlphaData() called");
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnModelAlphaData");

}

void HRLModelTestCommand::RunSimulationTestOnModelAlphaMember(void)
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::RunSimulationTestOnModelAlphaMember() called");
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnModelAlphaMember");

}
