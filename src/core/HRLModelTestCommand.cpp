#include "HRLModelTestCommand.hpp"
#include "HRLModelHelper.hpp"
#include "HRLModelTester.hpp"
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
        {"0", TesterType::BENCHMARK},         {"benchmark",      TesterType::BENCHMARK},
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
        case TesterType::BENCHMARK:
            RunSimulationTestOnBenchMark();
            break;
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
                        "  [0] Simulation Test on Benchmark\n"
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

void HRLModelTestCommand::RunSimulationTestOnBenchMark(void)
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::RunSimulationTestOnBenchMark() called");
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnBenchMark");

    const int gaus_par_size{ 3 };
    const int linear_basis_size{ 2 };
    Eigen::VectorXd model_par_prior{ Eigen::VectorXd::Zero(gaus_par_size) };
    model_par_prior(0) = 1.0;
    model_par_prior(1) = 0.5;
    model_par_prior(2) = 0.1;
    Eigen::VectorXd model_par_sigma{ Eigen::VectorXd::Zero(gaus_par_size) };
    model_par_sigma(0) = 0.05;
    model_par_sigma(1) = 0.025;
    model_par_sigma(2) = 0.001;

    Eigen::VectorXd residual_mean_ols;
    Eigen::VectorXd residual_mean_mdpde;
    Eigen::VectorXd residual_sigma_ols;
    Eigen::VectorXd residual_sigma_mdpde;

    auto replica_size{ 1000 };
    auto tester{ std::make_unique<HRLModelTester>(gaus_par_size, linear_basis_size, replica_size) };
    tester->SetFittingRange(m_options.fit_range_min, m_options.fit_range_max);
    tester->SetSamplingEntrySize(1000);
    tester->SetDataErrorSigma(0.1);
    tester->SetOutlierRatio(0.2);
    tester->SetGausParametersPrior(model_par_prior);
    tester->SetGausParametersSigma(model_par_sigma);
    tester->RunBetaMDPDETest(
        m_options.alpha_r,
        residual_mean_ols, residual_mean_mdpde, residual_sigma_ols, residual_sigma_mdpde,
        m_options.thread_size
    );

    std::cout <<"OLS:   "
              << residual_mean_ols(0) <<" +- " << residual_sigma_ols(0) <<" , "
              << residual_mean_ols(1) <<" +- " << residual_sigma_ols(1) <<" , "
              << residual_mean_ols(2) <<" +- " << residual_sigma_ols(2) << std::endl;
    std::cout <<"MDODE: "
              << residual_mean_mdpde(0) <<" +- " << residual_sigma_mdpde(0) <<" , "
              << residual_mean_mdpde(1) <<" +- " << residual_sigma_mdpde(1) <<" , "
              << residual_mean_mdpde(2) <<" +- " << residual_sigma_mdpde(2) << std::endl;
}

void HRLModelTestCommand::RunSimulationTestOnDataOutlier(void)
{
    Logger::Log(LogLevel::Debug, "HRLModelTestCommand::RunSimulationTestOnDataOutlier() called");
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnDataOutlier");

    const int gaus_par_size{ 3 };
    const int linear_basis_size{ 2 };
    Eigen::VectorXd model_par_prior{ Eigen::VectorXd::Zero(gaus_par_size) };
    model_par_prior(0) = 1.0;
    model_par_prior(1) = 0.5;
    model_par_prior(2) = 0.1;
    Eigen::VectorXd model_par_sigma{ Eigen::VectorXd::Zero(gaus_par_size) };
    model_par_sigma(0) = 0.05;
    model_par_sigma(1) = 0.025;
    model_par_sigma(2) = 0.001;

    auto replica_size{ 100 };
    auto tester{ std::make_unique<HRLModelTester>(gaus_par_size, linear_basis_size, replica_size) };
    tester->SetFittingRange(m_options.fit_range_min, m_options.fit_range_max);
    tester->SetSamplingEntrySize(1000);
    tester->SetGausParametersPrior(model_par_prior);
    tester->SetGausParametersSigma(model_par_sigma);

    std::vector<double> error_list{ 0.1, 0.2, 0.3 };
    std::vector<Eigen::MatrixXd> mean_matrix_ols_list;
    std::vector<Eigen::MatrixXd> mean_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_ols_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_mdpde_list;
    
    auto outlier_size{ 9 };
    auto outlier_step_size{ 0.025 };
    for (auto error_sigma : error_list)
    {
        Eigen::MatrixXd mean_matrix_ols{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd mean_matrix_mdpde{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_ols{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_mdpde{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        for (int i = 0; i < outlier_size; i++)
        {
            auto outlier_ratio{ i * outlier_step_size };
            Eigen::VectorXd residual_mean_ols;
            Eigen::VectorXd residual_mean_mdpde;
            Eigen::VectorXd residual_sigma_ols;
            Eigen::VectorXd residual_sigma_mdpde;
            tester->SetDataErrorSigma(error_sigma);
            tester->SetOutlierRatio(outlier_ratio);
            tester->RunBetaMDPDETest(
                m_options.alpha_r,
                residual_mean_ols, residual_mean_mdpde, residual_sigma_ols, residual_sigma_mdpde,
                m_options.thread_size
            );

            mean_matrix_ols.col(i) = residual_mean_ols;
            mean_matrix_mdpde.col(i) = residual_mean_mdpde;
            sigma_matrix_ols.col(i) = residual_sigma_ols;
            sigma_matrix_mdpde.col(i) = residual_sigma_mdpde;
        }
        mean_matrix_ols_list.emplace_back(mean_matrix_ols);
        mean_matrix_mdpde_list.emplace_back(mean_matrix_mdpde);
        sigma_matrix_ols_list.emplace_back(sigma_matrix_ols);
        sigma_matrix_mdpde_list.emplace_back(sigma_matrix_mdpde);
    }
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
