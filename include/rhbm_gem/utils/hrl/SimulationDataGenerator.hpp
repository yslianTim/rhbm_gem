#pragma once

#include <random>

#include <Eigen/Dense>

#include <rhbm_gem/utils/math/SamplingTypes.hpp>

class SimulationDataGenerator
{
public:
    struct NeighborhoodOptions
    {
        double radius_min{ 0.0 };
        double radius_max{ 1.0 };
        double neighbor_distance{ 2.0 };
        size_t neighbor_count{ 1 };
        double angle{ 0.0 };
    };

private:
    int m_gaus_par_size;
    double m_x_min, m_x_max;

public:
    SimulationDataGenerator() = delete;
    explicit SimulationDataGenerator(int gaus_par_size);
    ~SimulationDataGenerator() = default;

    void SetFittingRange(double x_min, double x_max) { m_x_min = x_min; m_x_max = x_max; }

    LocalPotentialSampleList GenerateGaussianSampling(
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        double outlier_ratio,
        std::mt19937 & generator
    ) const;
    LocalPotentialSampleList GenerateGaussianSamplingWithNeighborhood(
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        const NeighborhoodOptions & options
    ) const;
    SeriesPointList GenerateLinearDataset(
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        double error_sigma,
        double outlier_ratio,
        std::mt19937 & generator
    ) const;
    SeriesPointList GenerateLinearDatasetWithNeighborhood(
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        double error_sigma,
        const NeighborhoodOptions & options,
        std::mt19937 & generator
    ) const;

private:
    void ValidateGausParametersDimension(const Eigen::VectorXd & gaus_par) const;
};
