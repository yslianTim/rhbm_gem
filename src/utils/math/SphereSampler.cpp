#include <rhbm_gem/utils/math/SphereSampler.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <cmath>
#include <random>
#include <stdexcept>
#include <sstream>

SphereSampler::SphereSampler() :
    m_sampling_size{ 10 }, m_distance_min{ 0.0 }, m_distance_max{ 1.0 }
{

}

void SphereSampler::Print() const
{
    std::ostringstream oss;
    oss << "SphereSampler Configuration:\n"
        << " - Sampling size: " << m_sampling_size << '\n'
        << " - Distance range: ["
        << StringHelper::ToStringWithPrecision<double>(m_distance_min, 1) << ", "
        << StringHelper::ToStringWithPrecision<double>(m_distance_max, 1)
        << "] Angstrom.";
    Logger::Log(LogLevel::Info, oss.str());
}

std::vector<std::tuple<float, std::array<float, 3>>> SphereSampler::GenerateSamplingPoints(
    const std::array<float, 3> & reference_position,
    const std::array<float, 3> & axis_vector) const
{
    (void)axis_vector; // axis_vector is not used in SphereSampler

    std::vector<std::tuple<float, std::array<float, 3>>> out;
    RunRadiusUniformRandomSamplingMethod(reference_position, out);

    return out;
}

void SphereSampler::RunVolumeUniformRandomSamplingMethod(
    const std::array<float, 3> & reference_position,
    std::vector<std::tuple<float, std::array<float, 3>>> & out) const
{
    if (m_distance_min > m_distance_max)
    {
        throw std::invalid_argument("SphereSampler: distance minimum greater than maximum");
    }
    if (m_distance_min < 0.0 || m_distance_max < 0.0)
    {
        throw std::invalid_argument("SphereSampler: distance range cannot be negative");
    }
    out.resize(m_sampling_size);
    
    static thread_local std::mt19937 engine{ std::random_device{}() };
    std::uniform_real_distribution<float> dist_unit(0.0f, 1.0f);
    std::uniform_real_distribution<float> dist_phi(0.0f, static_cast<float>(Constants::two_pi));
    std::uniform_real_distribution<float> dist_cos_theta(-1.0f, 1.0f);
    const float min_radius_cube{
        static_cast<float>(m_distance_min * m_distance_min * m_distance_min)
    };
    const float max_radius_cube{
        static_cast<float>(m_distance_max * m_distance_max * m_distance_max)
    };

    for (unsigned int i = 0; i < m_sampling_size; i++)
    {
        const float radius_unit{ dist_unit(engine) };
        const float radius{
            std::cbrt(min_radius_cube + radius_unit * (max_radius_cube - min_radius_cube))
        };
        float phi{ dist_phi(engine) };
        float cos_theta{ dist_cos_theta(engine) };
        float sin_theta{ std::sqrt(1.0f - cos_theta * cos_theta) };

        std::array<float, 3> direction_vector{
            sin_theta * std::cos(phi),
            sin_theta * std::sin(phi),
            cos_theta
        };

        std::array<float, 3> sampling_position{
            reference_position.at(0) + radius * direction_vector[0],
            reference_position.at(1) + radius * direction_vector[1],
            reference_position.at(2) + radius * direction_vector[2]
        };

        out[i] = std::make_tuple(radius, sampling_position);
    }
}

void SphereSampler::RunRadiusUniformRandomSamplingMethod(
    const std::array<float, 3> & reference_position,
    std::vector<std::tuple<float, std::array<float, 3>>> & out) const
{
    if (m_distance_min > m_distance_max)
    {
        throw std::invalid_argument("SphereSampler: distance minimum greater than maximum");
    }
    if (m_distance_min < 0.0 || m_distance_max < 0.0)
    {
        throw std::invalid_argument("SphereSampler: distance range cannot be negative");
    }
    out.resize(m_sampling_size);
    
    static thread_local std::mt19937 engine{ std::random_device{}() };
    std::uniform_real_distribution<float> dist_radius(
        static_cast<float>(m_distance_min), static_cast<float>(m_distance_max)
    );
    std::uniform_real_distribution<float> dist_phi(0.0f, static_cast<float>(Constants::two_pi));
    std::uniform_real_distribution<float> dist_cos_theta(-1.0f, 1.0f);

    for (unsigned int i = 0; i < m_sampling_size; i++)
    {
        const float radius{ dist_radius(engine) };
        float phi{ dist_phi(engine) };
        float cos_theta{ dist_cos_theta(engine) };
        float sin_theta{ std::sqrt(1.0f - cos_theta * cos_theta) };

        std::array<float, 3> direction_vector{
            sin_theta * std::cos(phi),
            sin_theta * std::sin(phi),
            cos_theta
        };

        std::array<float, 3> sampling_position{
            reference_position.at(0) + radius * direction_vector[0],
            reference_position.at(1) + radius * direction_vector[1],
            reference_position.at(2) + radius * direction_vector[2]
        };

        out[i] = std::make_tuple(radius, sampling_position);
    }
}
