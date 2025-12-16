#pragma once

#include <CLI/CLI.hpp>
#include <filesystem>
#include <string>
#include <vector>
#include <Eigen/Dense>

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
    void PrintDataOutlierResult(
        const std::string & name,
        const std::vector<double> & outlier_list,
        const std::vector<Eigen::MatrixXd> & mean_matrix_ols_list,
        const std::vector<Eigen::MatrixXd> & mean_matrix_mdpde_list,
        const std::vector<Eigen::MatrixXd> & sigma_matrix_ols_list,
        const std::vector<Eigen::MatrixXd> & sigma_matrix_mdpde_list
    );
    void PrintMemberOutlierResult(
        const std::string & name,
        const std::vector<double> & outlier_list,
        const std::vector<Eigen::MatrixXd> & mean_matrix_median_list,
        const std::vector<Eigen::MatrixXd> & mean_matrix_mdpde_list,
        const std::vector<Eigen::MatrixXd> & sigma_matrix_median_list,
        const std::vector<Eigen::MatrixXd> & sigma_matrix_mdpde_list
    );

};
