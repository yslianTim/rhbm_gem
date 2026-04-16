#include <rhbm_gem/utils/hrl/SimulationDataGenerator.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>
#include <rhbm_gem/utils/math/SamplingPointFilter.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

SimulationDataGenerator::SimulationDataGenerator(int gaus_par_size) :
    m_gaus_par_size{ gaus_par_size },
    m_x_min{ 0.0 }, m_x_max{ 1.0 }
{
    if (m_gaus_par_size <= 0)
    {
        throw std::invalid_argument("gaus_par_size must be positive value");
    }
}

LocalPotentialSampleList SimulationDataGenerator::GenerateGaussianSampling(
    size_t sampling_entry_size,
    const Eigen::VectorXd & gaus_par,
    double outlier_ratio,
    std::mt19937 & generator
) const
{
    ValidateGausParametersDimension(gaus_par);
    const auto amplitude{ gaus_par(0) };
    const auto width{ gaus_par(1) };
    const auto intersect{ 0.0 };
    std::uniform_real_distribution<> dist_distance(m_x_min, m_x_max);
    std::uniform_real_distribution<> dist_outlier(0.0, 1.0);
    LocalPotentialSampleList sampling_entry_list;
    sampling_entry_list.reserve(sampling_entry_size);
    for (size_t i = 0; i < sampling_entry_size; i++)
    {
        const auto r{ dist_distance(generator) };
        auto y{
            amplitude * GausLinearTransformHelper::GetGaussianResponseAtDistance(r, width) + intersect
        };
        const auto y_outlier{
            0.5 * amplitude * std::pow(Constants::two_pi * std::pow(width, 2), -1.5)
        };
        if (dist_outlier(generator) < outlier_ratio)
        {
            y = y_outlier;
        }
        sampling_entry_list.emplace_back(LocalPotentialSample{
            static_cast<float>(r),
            static_cast<float>(y)
        });
    }
    return sampling_entry_list;
}

LocalPotentialSampleList SimulationDataGenerator::GenerateGaussianSamplingWithNeighborhood(
    size_t sampling_entry_size,
    const Eigen::VectorXd & gaus_par,
    const NeighborhoodOptions & options
) const
{
    ValidateGausParametersDimension(gaus_par);
    constexpr size_t max_neighbor_count{ 4 };
    if (options.neighbor_count > max_neighbor_count)
    {
        throw std::invalid_argument("neighbor_count should be less than or equal to 4");
    }

    const auto amplitude{ gaus_par(0) };
    const auto width{ gaus_par(1) };
    const auto intersect{ 0.0 };
    Eigen::VectorXd atom_center{ Eigen::VectorXd::Zero(3) };
    Eigen::VectorXd neighbor_center[max_neighbor_count];
    for (size_t i = 0; i < max_neighbor_count; i++)
    {
        neighbor_center[i] = Eigen::VectorXd::Zero(3);
    }

    if (options.neighbor_count <= 2)
    {
        neighbor_center[0] << 1.0, 0.0, 0.0;
        neighbor_center[1] << -1.0, 0.0, 0.0;
    }

    if (options.neighbor_count == 3)
    {
        neighbor_center[0] << 1.0, 0.0, 0.0;
        neighbor_center[1] << -0.5, std::sqrt(3) / 2.0, 0.0;
        neighbor_center[2] << -0.5, -std::sqrt(3) / 2.0, 0.0;
    }

    if (options.neighbor_count == 4)
    {
        neighbor_center[0] << 0.0, 0.0, 1.0;
        neighbor_center[1] << 0.0, 2.0 * std::sqrt(2) / 3.0, -1.0 / 3.0;
        neighbor_center[2] << -std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0;
        neighbor_center[3] << std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0;
    }

    std::vector<Eigen::VectorXd> neighbor_center_list;
    neighbor_center_list.reserve(options.neighbor_count);
    for (size_t i = 0; i < options.neighbor_count; i++)
    {
        neighbor_center[i] *= options.neighbor_distance;
        neighbor_center_list.emplace_back(neighbor_center[i]);
    }

    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ options.radius_min, options.radius_max },
            0.1,
            static_cast<unsigned int>(sampling_entry_size),
            false
        )
    );
    const auto sampling_points{ sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }) };
    const auto filtered_sampling_points{
        rhbm_gem::SelectSamplingPoint(sampling_points, neighbor_center_list, options.angle)
    };
    LocalPotentialSampleList sampling_entry_list;
    sampling_entry_list.reserve(filtered_sampling_points.size());
    for (const auto & sampling_point : filtered_sampling_points)
    {
        Eigen::VectorXd point{ Eigen::VectorXd::Zero(3) };
        point(0) = sampling_point.position[0];
        point(1) = sampling_point.position[1];
        point(2) = sampling_point.position[2];
        const auto y{
            amplitude * GausLinearTransformHelper::GetGaussianPesponseAtPointWithNeighborhood(
                point,
                atom_center,
                neighbor_center_list,
                width
            ) + intersect
        };
        const auto weight{ 1.0f };
        sampling_entry_list.emplace_back(LocalPotentialSample{
            sampling_point.distance,
            static_cast<float>(y),
            weight,
            sampling_point.position
        });
    }
    return sampling_entry_list;
}

SeriesPointList SimulationDataGenerator::GenerateLinearDataset(
    size_t sampling_entry_size,
    const Eigen::VectorXd & gaus_par,
    double error_sigma,
    double outlier_ratio,
    std::mt19937 & generator
) const
{
    const auto sampling_entries{
        GenerateGaussianSampling(sampling_entry_size, gaus_par, outlier_ratio, generator)
    };
    auto linear_data_entry_list{
        GausLinearTransformHelper::MapValueTransform(sampling_entries, m_x_min, m_x_max)
    };
    const auto max_response{
        gaus_par(0) * GausLinearTransformHelper::GetGaussianResponseAtDistance(0.0, gaus_par(1))
    };
    std::normal_distribution<> dist_error(0.0, error_sigma * max_response);
    for (auto & data_entry : linear_data_entry_list)
    {
        data_entry.response += dist_error(generator);
    }
    return linear_data_entry_list;
}

SeriesPointList SimulationDataGenerator::GenerateLinearDatasetWithNeighborhood(
    size_t sampling_entry_size,
    const Eigen::VectorXd & gaus_par,
    double error_sigma,
    const NeighborhoodOptions & options,
    std::mt19937 & generator
) const
{
    const auto sampling_entries{
        GenerateGaussianSamplingWithNeighborhood(sampling_entry_size, gaus_par, options)
    };
    auto linear_data_entry_list{
        GausLinearTransformHelper::MapValueTransform(sampling_entries, m_x_min, m_x_max)
    };
    const auto max_response{
        gaus_par(0) * GausLinearTransformHelper::GetGaussianResponseAtDistance(0.0, gaus_par(1))
    };
    std::normal_distribution<> dist_error(0.0, error_sigma * max_response);
    for (auto & data_entry : linear_data_entry_list)
    {
        data_entry.response += dist_error(generator);
    }
    return linear_data_entry_list;
}

void SimulationDataGenerator::ValidateGausParametersDimension(const Eigen::VectorXd & gaus_par) const
{
    if (gaus_par.rows() != m_gaus_par_size)
    {
        throw std::invalid_argument(
            "model parameters size invalid, must be : " + std::to_string(m_gaus_par_size)
        );
    }
}
