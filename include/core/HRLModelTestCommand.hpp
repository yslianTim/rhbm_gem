#pragma once

#include <CLI/CLI.hpp>
#include <filesystem>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "CommandBase.hpp"
#include "OptionEnumClass.hpp"

namespace rhbm_gem {

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
    HRLModelTestCommand();
    ~HRLModelTestCommand() = default;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    CommandId GetCommandId() const override { return CommandId::ModelTest; }
    const CommandOptions & GetOptions() const override { return m_options; }
    CommandOptions & GetOptions() override { return m_options; }

    void SetTesterChoice(TesterType value);
    void SetFitRangeMinimum(double value);
    void SetFitRangeMaximum(double value);
    void SetAlphaR(double value);
    void SetAlphaG(double value);

private:
    bool ExecuteImpl() override;
    void RunSimulationTestOnBenchMark();
    void RunSimulationTestOnDataOutlier();
    void RunSimulationTestOnMemberOutlier();
    void RunSimulationTestOnModelAlphaData();
    void RunSimulationTestOnModelAlphaMember();
    void PrintDataOutlierResult(
        const std::string & name,
        const std::vector<double> & outlier_list,
        const std::vector<Eigen::MatrixXd> & mean_matrix_ols_list,
        const std::vector<Eigen::MatrixXd> & mean_matrix_mdpde_list,
        const std::vector<Eigen::MatrixXd> & mean_matrix_train_list,
        const std::vector<Eigen::MatrixXd> & sigma_matrix_ols_list,
        const std::vector<Eigen::MatrixXd> & sigma_matrix_mdpde_list,
        const std::vector<Eigen::MatrixXd> & sigma_matrix_train_list
    );
    void PrintMemberOutlierResult(
        const std::string & name,
        const std::vector<double> & outlier_list,
        const std::vector<Eigen::MatrixXd> & mean_matrix_median_list,
        const std::vector<Eigen::MatrixXd> & mean_matrix_mdpde_list,
        const std::vector<Eigen::MatrixXd> & mean_matrix_train_list,
        const std::vector<Eigen::MatrixXd> & sigma_matrix_median_list,
        const std::vector<Eigen::MatrixXd> & sigma_matrix_mdpde_list,
        const std::vector<Eigen::MatrixXd> & sigma_matrix_train_list
    );

};

} // namespace rhbm_gem
