#pragma once

#include <cstddef>
#include <functional>
#include <sstream>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

class HRLAlphaTrainer
{
    double m_alpha_min{ 0.0 };
    double m_alpha_max{ 0.0 };
    double m_alpha_step{ 0.0 };
    std::vector<double> m_alpha_grid;

public:
    using ProgressCallback = std::function<void(std::size_t completed, std::size_t total)>;

    struct AlphaRunOptions
    {
        RHBMExecutionOptions execution_options{};
        ProgressCallback progress_callback{};
    };

    struct AlphaTrainingOptions : AlphaRunOptions
    {
        std::size_t subset_size{ 0 };
    };

    struct AlphaTrainingResult
    {
        double best_alpha{ 0.0 };
        Eigen::VectorXd error_sum_list;
    };

    HRLAlphaTrainer(double alpha_min, double alpha_max, double alpha_step);

    const std::vector<double> & AlphaGrid() const { return m_alpha_grid; }
    std::ostringstream GetAlphaGridSummary() const;

    AlphaTrainingResult TrainAlphaR(
        const std::vector<RHBMMemberDataset> & dataset_list,
        const AlphaTrainingOptions & options
    ) const;

    AlphaTrainingResult TrainAlphaG(
        const std::vector<std::vector<RHBMBetaVector>> & beta_group_list,
        const AlphaTrainingOptions & options
    ) const;

    Eigen::MatrixXd StudyAlphaRBias(
        const std::vector<RHBMMemberDataset> & dataset_list,
        const AlphaRunOptions & options
    ) const;

    Eigen::MatrixXd StudyAlphaGBias(
        const std::vector<std::vector<RHBMBetaVector>> & beta_group_list,
        const AlphaRunOptions & options
    ) const;

private:
    std::vector<double> BuildAlphaGrid(double alpha_min, double alpha_max, double alpha_step);
};
