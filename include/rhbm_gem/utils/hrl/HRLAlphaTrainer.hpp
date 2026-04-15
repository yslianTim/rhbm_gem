#pragma once

#include <cstddef>
#include <functional>
#include <sstream>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>

class HRLAlphaTrainer
{
    double m_alpha_min{ 0.0 };
    double m_alpha_max{ 0.0 };
    double m_alpha_step{ 0.0 };
    std::vector<double> m_alpha_grid;

public:
    using ProgressCallback = std::function<void(std::size_t completed, std::size_t total)>;

    struct AlphaRTrainingOptions
    {
        std::size_t subset_size{ 5 };
        HRLExecutionOptions execution_options{};
        ProgressCallback progress_callback{};
    };

    struct AlphaGTrainingOptions
    {
        std::size_t subset_size{ 10 };
        HRLExecutionOptions execution_options{};
        ProgressCallback progress_callback{};
    };

    struct AlphaTrainingResult
    {
        double best_alpha{ 0.0 };
        Eigen::VectorXd error_sum_list;
    };

    struct AlphaBiasStudyOptions
    {
        HRLExecutionOptions execution_options{};
        ProgressCallback progress_callback{};
    };

    HRLAlphaTrainer(double alpha_min, double alpha_max, double alpha_step);

    const std::vector<double> & AlphaGrid() const { return m_alpha_grid; }
    std::ostringstream GetAlphaGridSummary() const;

    static Eigen::VectorXd EvaluateAlphaR(
        const HRLMemberDataset & dataset,
        std::size_t subset_size,
        const std::vector<double> & alpha_list,
        const HRLExecutionOptions & options = {}
    );

    static Eigen::VectorXd EvaluateAlphaG(
        const std::vector<Eigen::VectorXd> & beta_list,
        std::size_t subset_size,
        const std::vector<double> & alpha_list,
        const HRLExecutionOptions & options = {}
    );

    AlphaTrainingResult TrainAlphaR(
        const std::vector<HRLMemberDataset> & dataset_list) const;

    AlphaTrainingResult TrainAlphaR(
        const std::vector<HRLMemberDataset> & dataset_list,
        const AlphaRTrainingOptions & options
    ) const;

    AlphaTrainingResult TrainAlphaG(
        const std::vector<std::vector<Eigen::VectorXd>> & beta_group_list) const;

    AlphaTrainingResult TrainAlphaG(
        const std::vector<std::vector<Eigen::VectorXd>> & beta_group_list,
        const AlphaGTrainingOptions & options
    ) const;

    Eigen::MatrixXd StudyAlphaRBias(
        const std::vector<HRLMemberDataset> & dataset_list) const;

    Eigen::MatrixXd StudyAlphaRBias(
        const std::vector<HRLMemberDataset> & dataset_list,
        const AlphaBiasStudyOptions & options
    ) const;

    Eigen::MatrixXd StudyAlphaGBias(
        const std::vector<std::vector<Eigen::VectorXd>> & beta_group_list) const;

    Eigen::MatrixXd StudyAlphaGBias(
        const std::vector<std::vector<Eigen::VectorXd>> & beta_group_list,
        const AlphaBiasStudyOptions & options
    ) const;

private:
    std::vector<double> BuildAlphaGrid(double alpha_min, double alpha_max, double alpha_step);
};
