#include "HRLModelTestCommand.hpp"
#include "HRLModelTester.hpp"
#include "ArrayStats.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"

#include <memory>
#include <random>
#include <sstream>
#include <vector>

#ifdef HAVE_ROOT
#include "ROOTHelper.hpp"
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#endif

namespace rhbm_gem {

void HRLModelTestCommand::RunSimulationTestOnBenchMark()
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnBenchMark");

    const int gaus_par_size{ 3 };
    const int linear_basis_size{ 2 };
    Eigen::VectorXd model_par_prior{ Eigen::VectorXd::Zero(gaus_par_size) };
    model_par_prior(0) = 1.0;
    model_par_prior(1) = 0.5;
    model_par_prior(2) = 0.1;

    std::vector<Eigen::VectorXd> residual_mean_ols_list;
    std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
    std::vector<Eigen::VectorXd> residual_sigma_ols_list;
    std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;

    std::vector<double> alpha_r_list{ 0.1 };
    auto replica_size{ 1000 };
    auto sampling_entry_size{ 1000 };
    auto error_sigma{ 1.0 };
    auto outlier_ratio{ 0.0 };
    auto tester{ std::make_unique<HRLModelTester>(gaus_par_size, linear_basis_size, replica_size) };
    tester->SetFittingRange(m_options.fit_range_min, m_options.fit_range_max);
    tester->RunBetaMDPDETest(
        alpha_r_list,
        residual_mean_ols_list, residual_mean_mdpde_list,
        residual_sigma_ols_list, residual_sigma_mdpde_list,
        model_par_prior, sampling_entry_size, error_sigma, outlier_ratio, m_options.thread_size
    );

    std::ostringstream ols_stream;
    ols_stream << "OLS: "
               << residual_mean_ols_list.front()(0) << " +- "
               << residual_sigma_ols_list.front()(0) << " , "
               << residual_mean_ols_list.front()(1) << " +- "
               << residual_sigma_ols_list.front()(1) << " , "
               << residual_mean_ols_list.front()(2) << " +- "
               << residual_sigma_ols_list.front()(2);
    Logger::Log(LogLevel::Info, ols_stream.str());

    std::ostringstream mdpde_stream;
    mdpde_stream << "MDPDE: "
                 << residual_mean_mdpde_list.front()(0) << " +- "
                 << residual_sigma_mdpde_list.front()(0) << " , "
                 << residual_mean_mdpde_list.front()(1) << " +- "
                 << residual_sigma_mdpde_list.front()(1) << " , "
                 << residual_mean_mdpde_list.front()(2) << " +- "
                 << residual_sigma_mdpde_list.front()(2);
    Logger::Log(LogLevel::Info, mdpde_stream.str());
}

void HRLModelTestCommand::RunSimulationTestOnDataOutlier()
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnDataOutlier");

    const int gaus_par_size{ 3 };
    const int linear_basis_size{ 2 };
    Eigen::VectorXd model_par_prior{ Eigen::VectorXd::Zero(gaus_par_size) };
    model_par_prior(0) = 1.0;
    model_par_prior(1) = 0.5;
    model_par_prior(2) = 0.1;

    auto replica_size{ 100 };
    auto sampling_entry_size{ 1000 };
    auto tester{ std::make_unique<HRLModelTester>(gaus_par_size, linear_basis_size, replica_size) };
    tester->SetFittingRange(m_options.fit_range_min, m_options.fit_range_max);

    std::vector<double> alpha_r_list{ 0.1 };
    std::vector<double> error_list{ 0.1, 0.2, 0.3 };
    std::vector<Eigen::MatrixXd> mean_matrix_ols_list;
    std::vector<Eigen::MatrixXd> mean_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> mean_matrix_train_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_ols_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_train_list;
    
    auto outlier_size{ 9 };
    auto outlier_step_size{ 0.025 };
    std::vector<double> outlier_list(static_cast<size_t>(outlier_size));
    for (size_t i = 0; i < static_cast<size_t>(outlier_size); i++)
    {
        outlier_list[i] = static_cast<double>(i) * outlier_step_size;
    }

    for (auto error_sigma : error_list)
    {
        Eigen::MatrixXd mean_matrix_ols{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd mean_matrix_mdpde{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd mean_matrix_train{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_ols{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_mdpde{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_train{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        for (int i = 0; i < outlier_size; i++)
        {
            std::vector<Eigen::VectorXd> residual_mean_ols_list;
            std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
            std::vector<Eigen::VectorXd> residual_sigma_ols_list;
            std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;
            tester->RunBetaMDPDETest(
                alpha_r_list,
                residual_mean_ols_list, residual_mean_mdpde_list,
                residual_sigma_ols_list, residual_sigma_mdpde_list,
                model_par_prior, sampling_entry_size, error_sigma,
                outlier_list[static_cast<size_t>(i)], m_options.thread_size
            );

            mean_matrix_ols.col(i) = residual_mean_ols_list.front();
            mean_matrix_mdpde.col(i) = residual_mean_mdpde_list.front();
            sigma_matrix_ols.col(i) = residual_sigma_ols_list.front();
            sigma_matrix_mdpde.col(i) = residual_sigma_mdpde_list.front();
            mean_matrix_train.col(i) = residual_mean_mdpde_list.back();
            sigma_matrix_train.col(i) = residual_sigma_mdpde_list.back();
        }
        mean_matrix_ols_list.emplace_back(mean_matrix_ols);
        mean_matrix_mdpde_list.emplace_back(mean_matrix_mdpde);
        mean_matrix_train_list.emplace_back(mean_matrix_train);
        sigma_matrix_ols_list.emplace_back(sigma_matrix_ols);
        sigma_matrix_mdpde_list.emplace_back(sigma_matrix_mdpde);
        sigma_matrix_train_list.emplace_back(sigma_matrix_train);
    }

    PrintDataOutlierResult(
        "bias_outlier_in_data.pdf",
        outlier_list,
        mean_matrix_ols_list, mean_matrix_mdpde_list, mean_matrix_train_list,
        sigma_matrix_ols_list, sigma_matrix_mdpde_list, sigma_matrix_train_list
    );
}

void HRLModelTestCommand::RunSimulationTestOnMemberOutlier()
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnMemberOutlier");

    const int gaus_par_size{ 3 };
    const int linear_basis_size{ 2 };
    Eigen::VectorXd model_par_prior{ Eigen::VectorXd::Zero(gaus_par_size) };
    Eigen::VectorXd model_par_sigma{ Eigen::VectorXd::Zero(gaus_par_size) };
    model_par_prior(0) = 1.0;
    model_par_prior(1) = 0.5;
    model_par_prior(2) = 0.1;
    model_par_sigma(0) = 0.050;
    model_par_sigma(1) = 0.025;
    model_par_sigma(2) = 0.010;

    auto replica_size{ 100 };
    auto member_size{ 100 };
    auto tester{ std::make_unique<HRLModelTester>(gaus_par_size, linear_basis_size, replica_size) };
    tester->SetFittingRange(m_options.fit_range_min, m_options.fit_range_max);

    std::vector<double> alpha_g_list{ 0.2 };
    std::vector<Eigen::VectorXd> outlier_prior_list{
        Eigen::VectorXd::Zero(gaus_par_size),
        Eigen::VectorXd::Zero(gaus_par_size)
    };
    outlier_prior_list.at(0)(0) = 1.50;
    outlier_prior_list.at(0)(1) = 0.50;
    outlier_prior_list.at(0)(2) = 0.10;
    outlier_prior_list.at(1)(0) = 1.00;
    outlier_prior_list.at(1)(1) = 1.00;
    outlier_prior_list.at(1)(2) = 0.10;

    std::vector<Eigen::MatrixXd> mean_matrix_median_list;
    std::vector<Eigen::MatrixXd> mean_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> mean_matrix_train_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_median_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_train_list;
    
    auto outlier_size{ 9 };
    auto outlier_step_size{ 0.025 };
    std::vector<double> outlier_list(static_cast<size_t>(outlier_size));
    for (size_t i = 0; i < static_cast<size_t>(outlier_size); i++)
    {
        outlier_list[i] = static_cast<double>(i) * outlier_step_size;
    }

    for (auto outlier_prior : outlier_prior_list)
    {
        Eigen::MatrixXd mean_matrix_median{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd mean_matrix_mdpde{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd mean_matrix_train{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_median{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_mdpde{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_train{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        for (int i = 0; i < outlier_size; i++)
        {
            std::vector<Eigen::VectorXd> residual_mean_median_list;
            std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
            std::vector<Eigen::VectorXd> residual_sigma_median_list;
            std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;
            tester->RunMuMDPDETest(
                alpha_g_list,
                residual_mean_median_list, residual_mean_mdpde_list,
                residual_sigma_median_list, residual_sigma_mdpde_list,
                member_size, model_par_prior, model_par_sigma,
                outlier_prior, model_par_sigma,
                outlier_list[static_cast<size_t>(i)], m_options.thread_size
            );

            mean_matrix_median.col(i) = residual_mean_median_list.front();
            mean_matrix_mdpde.col(i) = residual_mean_mdpde_list.front();
            mean_matrix_train.col(i) = residual_mean_mdpde_list.back();
            sigma_matrix_median.col(i) = residual_sigma_median_list.front();
            sigma_matrix_mdpde.col(i) = residual_sigma_mdpde_list.front();
            sigma_matrix_train.col(i) = residual_sigma_mdpde_list.back();
        }
        mean_matrix_median_list.emplace_back(mean_matrix_median);
        mean_matrix_mdpde_list.emplace_back(mean_matrix_mdpde);
        mean_matrix_train_list.emplace_back(mean_matrix_train);
        sigma_matrix_median_list.emplace_back(sigma_matrix_median);
        sigma_matrix_mdpde_list.emplace_back(sigma_matrix_mdpde);
        sigma_matrix_train_list.emplace_back(sigma_matrix_train);
    }

    PrintMemberOutlierResult(
        "bias_outlier_in_member.pdf",
        outlier_list,
        mean_matrix_median_list, mean_matrix_mdpde_list, mean_matrix_train_list,
        sigma_matrix_median_list, sigma_matrix_mdpde_list, sigma_matrix_train_list
    );
}

void HRLModelTestCommand::RunSimulationTestOnModelAlphaData()
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnModelAlphaData");

    const int gaus_par_size{ 3 };
    const int linear_basis_size{ 2 };
    Eigen::VectorXd model_par_prior{ Eigen::VectorXd::Zero(gaus_par_size) };
    model_par_prior(0) = 1.0;
    model_par_prior(1) = 0.5;
    model_par_prior(2) = 0.1;

    auto replica_size{ 100 };
    auto sampling_entry_size{ 1000 };
    auto tester{ std::make_unique<HRLModelTester>(gaus_par_size, linear_basis_size, replica_size) };
    tester->SetFittingRange(m_options.fit_range_min, m_options.fit_range_max);

    std::vector<double> alpha_r_list{ 0.1, 0.4 };
    std::vector<double> error_list{ 0.1, 0.2, 0.3 };
    std::vector<Eigen::MatrixXd> mean_matrix_alpha1_list;
    std::vector<Eigen::MatrixXd> mean_matrix_alpha2_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_alpha1_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_alpha2_list;
    
    auto outlier_size{ 10 };
    auto outlier_step_size{ 0.05 };
    std::vector<double> outlier_list(static_cast<size_t>(outlier_size));
    for (size_t i = 0; i < static_cast<size_t>(outlier_size); i++)
    {
        outlier_list[i] = static_cast<double>(i) * outlier_step_size;
    }

    for (auto error_sigma : error_list)
    {
        Eigen::MatrixXd mean_matrix_alpha1{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd mean_matrix_alpha2{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_alpha1{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_alpha2{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        for (int i = 0; i < outlier_size; i++)
        {
            std::vector<Eigen::VectorXd> residual_mean_ols_list;
            std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
            std::vector<Eigen::VectorXd> residual_sigma_ols_list;
            std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;
            tester->RunBetaMDPDETest(
                alpha_r_list,
                residual_mean_ols_list, residual_mean_mdpde_list,
                residual_sigma_ols_list, residual_sigma_mdpde_list,
                model_par_prior, sampling_entry_size, error_sigma,
                outlier_list[static_cast<size_t>(i)], m_options.thread_size
            );

            mean_matrix_alpha1.col(i) = residual_mean_mdpde_list.front();
            mean_matrix_alpha2.col(i) = residual_mean_mdpde_list.back();
            sigma_matrix_alpha1.col(i) = residual_sigma_mdpde_list.front();
            sigma_matrix_alpha2.col(i) = residual_sigma_mdpde_list.back();
        }
        mean_matrix_alpha1_list.emplace_back(mean_matrix_alpha1);
        mean_matrix_alpha2_list.emplace_back(mean_matrix_alpha2);
        sigma_matrix_alpha1_list.emplace_back(sigma_matrix_alpha1);
        sigma_matrix_alpha2_list.emplace_back(sigma_matrix_alpha2);
    }

    PrintDataOutlierResult(
        "bias_outlier_with_alpha_in_data.pdf",
        outlier_list,
        mean_matrix_alpha1_list, mean_matrix_alpha2_list, mean_matrix_alpha2_list,
        sigma_matrix_alpha1_list, sigma_matrix_alpha2_list, sigma_matrix_alpha2_list
    );
}

void HRLModelTestCommand::RunSimulationTestOnModelAlphaMember()
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnModelAlphaMember");

    const int gaus_par_size{ 3 };
    const int linear_basis_size{ 2 };
    Eigen::VectorXd model_par_prior{ Eigen::VectorXd::Zero(gaus_par_size) };
    Eigen::VectorXd model_par_sigma{ Eigen::VectorXd::Zero(gaus_par_size) };
    model_par_prior(0) = 1.0;
    model_par_prior(1) = 0.5;
    model_par_prior(2) = 0.1;
    model_par_sigma(0) = 0.050;
    model_par_sigma(1) = 0.025;
    model_par_sigma(2) = 0.010;

    auto replica_size{ 100 };
    auto member_size{ 100 };
    auto tester{ std::make_unique<HRLModelTester>(gaus_par_size, linear_basis_size, replica_size) };
    tester->SetFittingRange(m_options.fit_range_min, m_options.fit_range_max);

    std::vector<double> alpha_g_list{ 0.2, 0.5 };
    std::vector<Eigen::VectorXd> outlier_prior_list{
        Eigen::VectorXd::Zero(gaus_par_size),
        Eigen::VectorXd::Zero(gaus_par_size)
    };
    outlier_prior_list.at(0)(0) = 1.50;
    outlier_prior_list.at(0)(1) = 0.50;
    outlier_prior_list.at(0)(2) = 0.10;
    outlier_prior_list.at(1)(0) = 1.00;
    outlier_prior_list.at(1)(1) = 1.00;
    outlier_prior_list.at(1)(2) = 0.10;

    std::vector<Eigen::MatrixXd> mean_matrix_alpha1_list;
    std::vector<Eigen::MatrixXd> mean_matrix_alpha2_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_alpha1_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_alpha2_list;
    
    auto outlier_size{ 10 };
    auto outlier_step_size{ 0.05 };
    std::vector<double> outlier_list(static_cast<size_t>(outlier_size));
    for (size_t i = 0; i < static_cast<size_t>(outlier_size); i++)
    {
        outlier_list[i] = static_cast<double>(i) * outlier_step_size;
    }

    for (auto outlier_prior : outlier_prior_list)
    {
        Eigen::MatrixXd mean_matrix_alpha1{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd mean_matrix_alpha2{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_alpha1{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        Eigen::MatrixXd sigma_matrix_alpha2{ Eigen::MatrixXd::Zero(gaus_par_size, outlier_size) };
        for (int i = 0; i < outlier_size; i++)
        {
            std::vector<Eigen::VectorXd> residual_mean_median_list;
            std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
            std::vector<Eigen::VectorXd> residual_sigma_median_list;
            std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;
            tester->RunMuMDPDETest(
                alpha_g_list,
                residual_mean_median_list, residual_mean_mdpde_list,
                residual_sigma_median_list, residual_sigma_mdpde_list,
                member_size, model_par_prior, model_par_sigma,
                outlier_prior, model_par_sigma,
                outlier_list[static_cast<size_t>(i)], m_options.thread_size
            );

            mean_matrix_alpha1.col(i) = residual_mean_mdpde_list.front();
            mean_matrix_alpha2.col(i) = residual_mean_mdpde_list.back();
            sigma_matrix_alpha1.col(i) = residual_sigma_mdpde_list.front();
            sigma_matrix_alpha2.col(i) = residual_sigma_mdpde_list.back();
        }
        mean_matrix_alpha1_list.emplace_back(mean_matrix_alpha1);
        mean_matrix_alpha2_list.emplace_back(mean_matrix_alpha2);
        sigma_matrix_alpha1_list.emplace_back(sigma_matrix_alpha1);
        sigma_matrix_alpha2_list.emplace_back(sigma_matrix_alpha2);
    }

    PrintMemberOutlierResult(
        "bias_outlier_with_alpha_in_member.pdf",
        outlier_list,
        mean_matrix_alpha1_list, mean_matrix_alpha2_list, mean_matrix_alpha2_list,
        sigma_matrix_alpha1_list, sigma_matrix_alpha2_list, sigma_matrix_alpha2_list
    );
}

void HRLModelTestCommand::PrintDataOutlierResult(
    const std::string & name,
    const std::vector<double> & outlier_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_ols_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_train_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_ols_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_train_list)
{
    auto file_path{ BuildOutputPath(name, "") };
    Logger::Log(LogLevel::Info, " HRLModelTestCommand::PrintDataOutlierResult");

    std::vector<std::string> title_y_list{
        "Amplitude #font[2]{A}", "Width #tau"
    };

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 3 };
    const int row_size{ 2 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 750) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.17f, 0.08f, 0.12f, 0.14f, 0.01f, 0.01f
    );
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::vector<std::unique_ptr<TGraphErrors>> graph_ols_list[col_size][row_size];
    std::vector<std::unique_ptr<TGraphErrors>> graph_mdpde_list[col_size][row_size];
    std::vector<std::unique_ptr<TGraphErrors>> graph_train_list[col_size][row_size];
    std::vector<double> y_array[col_size][row_size];
    std::vector<double> global_y_array[row_size];
    for (size_t i = 0; i < col_size; i++)
    {
        for (size_t j = 0; j < row_size; j++)
        {
            auto graph_ols{ ROOTHelper::CreateGraphErrors() };
            auto graph_mdpde{ ROOTHelper::CreateGraphErrors() };
            auto graph_train{ ROOTHelper::CreateGraphErrors() };
            for (int p = 0; p < static_cast<int>(outlier_list.size()); p++)
            {
                auto x_value{ outlier_list.at(static_cast<size_t>(p)) * 100.0 };
                auto mean_ols{ mean_matrix_ols_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_ols{ sigma_matrix_ols_list.at(i).col(p)(static_cast<int>(j)) };
                auto mean_mdpde{ mean_matrix_mdpde_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_mdpde{ sigma_matrix_mdpde_list.at(i).col(p)(static_cast<int>(j)) };
                auto mean_train{ mean_matrix_train_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_train{ sigma_matrix_train_list.at(i).col(p)(static_cast<int>(j)) };
                graph_ols->SetPoint(p, x_value, mean_ols);
                graph_ols->SetPointError(p, 0.0, sigma_ols);
                graph_mdpde->SetPoint(p, x_value, mean_mdpde);
                graph_mdpde->SetPointError(p, 0.0, sigma_mdpde);
                graph_train->SetPoint(p, x_value, mean_train);
                graph_train->SetPointError(p, 0.0, sigma_train);
                y_array[i][j].emplace_back(mean_ols);
                y_array[i][j].emplace_back(mean_mdpde);
                y_array[i][j].emplace_back(mean_train);
                global_y_array[i].emplace_back(mean_ols);
                global_y_array[j].emplace_back(mean_mdpde);
                global_y_array[j].emplace_back(mean_train);
            }
            graph_ols_list[i][j].emplace_back(std::move(graph_ols));
            graph_mdpde_list[i][j].emplace_back(std::move(graph_mdpde));
            graph_train_list[i][j].emplace_back(std::move(graph_train));
        }
    }

    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 0.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 0.0 };
    for (size_t i = 0; i < col_size; i++)
    {
        x_min[i] = (m_options.tester_choice == TesterType::MODEL_ALPHA_DATA) ? -2.0 : -0.7;
        x_max[i] = (m_options.tester_choice == TesterType::MODEL_ALPHA_DATA) ? 47.0 : 22.0;
    }
    for (size_t j = 0; j < row_size; j++)
    {
        auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(global_y_array[j], 0.38, 0.005, 0.995) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> resolution_text[col_size];
    std::unique_ptr<TPaveText> title_x_text[col_size];
    std::unique_ptr<TPaveText> title_y_text[row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto par_id{ row_size - j - 1 };
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min[par_id], y_max[par_id]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 50.0f, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.08/x_factor), 505);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 50.0f, 1.2f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw("Y+");

            short color_ols{ kAzure };
            short color_mdpde{ kRed };
            short color_train{ kGreen+1 };
            for (auto & graph : graph_ols_list[i][par_id])
            {
                ROOTHelper::SetMarkerAttribute(graph.get(), 24, 1.5f, color_ols);
                ROOTHelper::SetLineAttribute(graph.get(), 2, 2, color_ols);
                ROOTHelper::SetFillAttribute(graph.get(), 1001, color_ols, 0.2f);
                graph->Draw("PL3");
            }
            for (auto & graph : graph_mdpde_list[i][par_id])
            {
                ROOTHelper::SetMarkerAttribute(graph.get(), 20, 1.5f, color_mdpde);
                ROOTHelper::SetLineAttribute(graph.get(), 1, 2, color_mdpde);
                ROOTHelper::SetFillAttribute(graph.get(), 1001, color_mdpde, 0.2f);
                graph->Draw("PL3");
            }
            for (auto & graph : graph_train_list[i][par_id])
            {
                ROOTHelper::SetMarkerAttribute(graph.get(), 25, 1.5f, color_train);
                ROOTHelper::SetLineAttribute(graph.get(), 3, 2, color_train);
                ROOTHelper::SetFillAttribute(graph.get(), 1001, color_train, 0.2f);
                graph->Draw("PL3");
            }

            if (i == 0)
            {
                title_y_text[j] = ROOTHelper::CreatePaveText(-0.68, 0.30, 0.00, 0.70, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(title_y_text[j].get());
                ROOTHelper::SetPaveAttribute(title_y_text[j].get(), 0, 0.1);
                ROOTHelper::SetLineAttribute(title_y_text[j].get(), 1, 0);
                ROOTHelper::SetTextAttribute(title_y_text[j].get(), 45.0f, 133, 22);
                ROOTHelper::SetFillAttribute(title_y_text[j].get(), 1001, kAzure-7, 0.5f);
                title_y_text[j]->AddText(title_y_list[static_cast<size_t>(par_id)].data());
                title_y_text[j]->Draw();
            }

            if (j == row_size - 1)
            {
                title_x_text[i] = ROOTHelper::CreatePaveText(0.10, 0.95, 0.90, 1.15, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(title_x_text[i].get());
                ROOTHelper::SetPaveAttribute(title_x_text[i].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(title_x_text[i].get(), 45.0f, 133, 22);
                ROOTHelper::SetFillAttribute(title_x_text[i].get(), 1001, kRed+1, 0.5f);
                title_x_text[i]->AddText(Form("#sigma_{#epsilon} = %.2f D_{max}", (i+1) * 0.1));
                title_x_text[i]->Draw();
            }
        }
    }

    canvas->cd();
    auto pad_extra0{ ROOTHelper::CreatePad("pad_extra0","", 0.02, 0.92, 0.98, 1.00) };
    pad_extra0->Draw();
    pad_extra0->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra0.get());
    ROOTHelper::SetFillAttribute(pad_extra0.get(), 4000);
    auto legend{ ROOTHelper::CreateLegend(0.0, 0.0, 1.0, 1.0, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12, 0.0);
    legend->SetMargin(0.25f);
    legend->SetNColumns(3);
    if (m_options.tester_choice == TesterType::MODEL_ALPHA_DATA)
    {
        legend->AddEntry(graph_ols_list[0][0].front().get(),
            "MDPDE (#alpha_{r} = 0.1)", "plf");
        legend->AddEntry(graph_mdpde_list[0][0].front().get(),
            "MDPDE (#alpha_{r} = 0.4)", "plf");
    }
    else
    {
        legend->AddEntry(graph_mdpde_list[0][0].front().get(),
            "MDPDE (#alpha_{r} = 0.1)", "plf");
        legend->AddEntry(graph_train_list[0][0].front().get(),
            "MDPDE (#alpha_{r} from training)", "plf");
        legend->AddEntry(graph_ols_list[0][0].front().get(),
            "Ordinary Least Squares", "plf");
    }
    legend->Draw();

    canvas->cd();
    auto pad_extra1{ ROOTHelper::CreatePad("pad_extra1","", 0.17, 0.00, 0.92, 0.06) };
    pad_extra1->Draw();
    pad_extra1->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra1.get());
    ROOTHelper::SetFillAttribute(pad_extra1.get(), 4000);
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 45.0f, 133, 22);
    bottom_title_text->AddText("Data Outlier Ratio (%)");
    bottom_title_text->Draw();

    canvas->cd();
    auto pad_extra2{ ROOTHelper::CreatePad("pad_extra2","", 0.96, 0.10, 1.00, 0.86) };
    pad_extra2->Draw();
    pad_extra2->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra2.get());
    ROOTHelper::SetFillAttribute(pad_extra2.get(), 4000);
    auto right_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(right_title_text.get());
    ROOTHelper::SetFillAttribute(right_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(right_title_text.get(), 50.0f, 133, 22);
    right_title_text->AddText("Normalized Bias #font[2]{b}");
    auto text{ right_title_text->GetLineWith("Bias") };
    text->SetTextAngle(90.0f);
    right_title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path.string());
    #endif
}

void HRLModelTestCommand::PrintMemberOutlierResult(
    const std::string & name,
    const std::vector<double> & outlier_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_ols_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_train_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_ols_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_train_list)
{
    auto file_path{ BuildOutputPath(name, "") };
    Logger::Log(LogLevel::Info, " HRLModelTestCommand::PrintMemberOutlierResult");

    std::vector<std::string> title_y_list{
        "Amplitude #font[2]{A}", "Width #tau"
    };
    std::vector<std::string> outlier_type_list{
        "#font[2]{A}", "#tau"
    };

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 2 };
    const int row_size{ 2 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 750) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.19f, 0.10f, 0.12f, 0.14f, 0.01f, 0.01f
    );
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::vector<std::unique_ptr<TGraphErrors>> graph_ols_list[col_size][row_size];
    std::vector<std::unique_ptr<TGraphErrors>> graph_mdpde_list[col_size][row_size];
    std::vector<std::unique_ptr<TGraphErrors>> graph_train_list[col_size][row_size];
    std::vector<double> global_y_array;
    for (size_t i = 0; i < col_size; i++)
    {
        for (size_t j = 0; j < row_size; j++)
        {
            auto graph_ols{ ROOTHelper::CreateGraphErrors() };
            auto graph_mdpde{ ROOTHelper::CreateGraphErrors() };
            auto graph_train{ ROOTHelper::CreateGraphErrors() };
            for (int p = 0; p < static_cast<int>(outlier_list.size()); p++)
            {
                auto x_value{ outlier_list.at(static_cast<size_t>(p)) * 100.0 };
                auto mean_ols{ mean_matrix_ols_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_ols{ sigma_matrix_ols_list.at(i).col(p)(static_cast<int>(j)) };
                auto mean_mdpde{ mean_matrix_mdpde_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_mdpde{ sigma_matrix_mdpde_list.at(i).col(p)(static_cast<int>(j)) };
                auto mean_train{ mean_matrix_train_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_train{ sigma_matrix_train_list.at(i).col(p)(static_cast<int>(j)) };
                graph_ols->SetPoint(p, x_value, mean_ols);
                graph_ols->SetPointError(p, 0.0, sigma_ols);
                graph_mdpde->SetPoint(p, x_value, mean_mdpde);
                graph_mdpde->SetPointError(p, 0.0, sigma_mdpde);
                graph_train->SetPoint(p, x_value, mean_train);
                graph_train->SetPointError(p, 0.0, sigma_train);
                global_y_array.emplace_back(mean_ols);
                global_y_array.emplace_back(mean_mdpde);
                global_y_array.emplace_back(mean_train);
            }
            graph_ols_list[i][j].emplace_back(std::move(graph_ols));
            graph_mdpde_list[i][j].emplace_back(std::move(graph_mdpde));
            graph_train_list[i][j].emplace_back(std::move(graph_train));
        }
    }

    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 0.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 0.0 };
    for (size_t i = 0; i < col_size; i++)
    {
        x_min[i] = (m_options.tester_choice == TesterType::MODEL_ALPHA_MEMBER) ? -2.0 : -0.7;
        x_max[i] = (m_options.tester_choice == TesterType::MODEL_ALPHA_MEMBER) ? 47.0 : 22.0;
    }
    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(global_y_array, 0.3) };
    for (size_t j = 0; j < row_size; j++)
    {
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> resolution_text[col_size];
    std::unique_ptr<TPaveText> title_x_text[col_size];
    std::unique_ptr<TPaveText> title_y_text[row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto par_id{ row_size - j - 1 };
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min[par_id], y_max[par_id]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 50.0f, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.08/x_factor), 505);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 50.0f, 1.2f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw("Y+");

            short color_ols{ kAzure };
            short color_mdpde{ kRed };
            short color_train{ kGreen+1 };
            for (auto & graph : graph_ols_list[i][par_id])
            {
                ROOTHelper::SetMarkerAttribute(graph.get(), 24, 1.5f, color_ols);
                ROOTHelper::SetLineAttribute(graph.get(), 2, 2, color_ols);
                ROOTHelper::SetFillAttribute(graph.get(), 1001, color_ols, 0.2f);
                graph->Draw("PL3");
            }
            for (auto & graph : graph_mdpde_list[i][par_id])
            {
                ROOTHelper::SetMarkerAttribute(graph.get(), 20, 1.5f, color_mdpde);
                ROOTHelper::SetLineAttribute(graph.get(), 1, 2, color_mdpde);
                ROOTHelper::SetFillAttribute(graph.get(), 1001, color_mdpde, 0.2f);
                graph->Draw("PL3");
            }
            for (auto & graph : graph_train_list[i][par_id])
            {
                ROOTHelper::SetMarkerAttribute(graph.get(), 25, 1.5f, color_train);
                ROOTHelper::SetLineAttribute(graph.get(), 3, 2, color_train);
                ROOTHelper::SetFillAttribute(graph.get(), 1001, color_train, 0.2f);
                graph->Draw("PL3");
            }

            if (i == 0)
            {
                title_y_text[j] = ROOTHelper::CreatePaveText(-0.52, 0.30, 0.00, 0.70, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(title_y_text[j].get());
                ROOTHelper::SetPaveAttribute(title_y_text[j].get(), 0, 0.1);
                ROOTHelper::SetLineAttribute(title_y_text[j].get(), 1, 0);
                ROOTHelper::SetTextAttribute(title_y_text[j].get(), 45.0f, 133, 22);
                ROOTHelper::SetFillAttribute(title_y_text[j].get(), 1001, kAzure-7, 0.5f);
                title_y_text[j]->AddText(title_y_list[static_cast<size_t>(par_id)].data());
                title_y_text[j]->Draw();
            }

            if (j == row_size - 1)
            {
                title_x_text[i] = ROOTHelper::CreatePaveText(0.10, 0.95, 0.90, 1.15, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(title_x_text[i].get());
                ROOTHelper::SetPaveAttribute(title_x_text[i].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(title_x_text[i].get(), 45.0f, 133, 22);
                ROOTHelper::SetFillAttribute(title_x_text[i].get(), 1001, kRed+1, 0.5f);
                title_x_text[i]->AddText(
                    Form("Outlier in %s", outlier_type_list[static_cast<size_t>(i)].data())
                );
                title_x_text[i]->Draw();
            }
        }
    }

    canvas->cd();
    auto pad_extra0{ ROOTHelper::CreatePad("pad_extra0","", 0.02, 0.92, 0.98, 1.00) };
    pad_extra0->Draw();
    pad_extra0->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra0.get());
    ROOTHelper::SetFillAttribute(pad_extra0.get(), 4000);
    auto legend{ ROOTHelper::CreateLegend(0.0, 0.0, 1.0, 1.0, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12, 0.0);
    legend->SetMargin(0.25f);
    legend->SetNColumns(3);
    if (m_options.tester_choice == TesterType::MODEL_ALPHA_MEMBER)
    {
        legend->AddEntry(graph_ols_list[0][0].front().get(),
            "MDPDE (#alpha_{g} = 0.2)", "plf");
        legend->AddEntry(graph_mdpde_list[0][0].front().get(),
            "MDPDE (#alpha_{g} = 0.5)", "plf");
    }
    else
    {
        legend->AddEntry(graph_mdpde_list[0][0].front().get(),
            "MDPDE (#alpha_{g} = 0.2)", "plf");
        legend->AddEntry(graph_train_list[0][0].front().get(),
            "MDPDE (#alpha_{g} from training)", "plf");
        legend->AddEntry(graph_ols_list[0][0].front().get(),
            "MDPDE (#alpha_{g} = 0)", "plf");
    }
    legend->Draw();

    canvas->cd();
    auto pad_extra1{ ROOTHelper::CreatePad("pad_extra1","", 0.19, 0.00, 0.90, 0.06) };
    pad_extra1->Draw();
    pad_extra1->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra1.get());
    ROOTHelper::SetFillAttribute(pad_extra1.get(), 4000);
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 45.0f, 133, 22);
    bottom_title_text->AddText("Member Outlier Ratio (%)");
    bottom_title_text->Draw();

    canvas->cd();
    auto pad_extra2{ ROOTHelper::CreatePad("pad_extra2","", 0.96, 0.10, 1.00, 0.86) };
    pad_extra2->Draw();
    pad_extra2->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra2.get());
    ROOTHelper::SetFillAttribute(pad_extra2.get(), 4000);
    auto right_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(right_title_text.get());
    ROOTHelper::SetFillAttribute(right_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(right_title_text.get(), 50.0f, 133, 22);
    right_title_text->AddText("Normalized Bias #font[2]{b}");
    auto text{ right_title_text->GetLineWith("Bias") };
    text->SetTextAngle(90.0f);
    right_title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path.string());
    #endif
}


} // namespace rhbm_gem
