#include "SphereSampler.hpp"
#include "Constants.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

#include <random>

#ifdef USE_OPENMP
#include <omp.h>
#endif

SphereSampler::SphereSampler(void) :
    m_thread_size{ 1 }, m_sampling_size{ 10 }, m_distance_min{ 0.0 }, m_distance_max{ 1.0 }
{

}

void SphereSampler::Print(void) const
{
    Logger::Log(LogLevel::Info,
                "SphereSampler Configuration:\n"
                " - Thread size: " + std::to_string(m_thread_size) + "\n"
                " - Sampling size: " + std::to_string(m_sampling_size) + "\n"
                " - Distance range: ["
                + StringHelper::ToStringWithPrecision<double>(m_distance_min, 1) + ", "
                + StringHelper::ToStringWithPrecision<double>(m_distance_max, 1) + "] Angstrom.");
}

std::vector<std::tuple<float, std::array<float, 3>>> SphereSampler::GenerateSamplingPoints(
    const std::array<float, 3> & position) const
{
    std::vector<std::tuple<float, std::array<float, 3>>> sampling_position_list;
    sampling_position_list.resize(m_sampling_size);

    #ifdef USE_OPENMP
    #pragma omp parallel for num_threads(m_thread_size)
    #endif
    for (unsigned int i = 0; i < m_sampling_size; i++)
    {
        static thread_local std::mt19937 engine{ std::random_device{}() };
        std::uniform_real_distribution<float> dist_r(static_cast<float>(m_distance_min), static_cast<float>(m_distance_max));
        std::uniform_real_distribution<float> dist_phi(0.0f, static_cast<float>(Constants::two_pi));
        std::uniform_real_distribution<float> dist_cos(-1.0f, 1.0f);

        float r{ dist_r(engine) };
        float phi{ dist_phi(engine) };
        float cos_theta{ dist_cos(engine) };
        float sin_theta{ std::sqrt(1.0f - cos_theta * cos_theta) };

        std::array<float, 3> direction{
            sin_theta * std::cos(phi),
            sin_theta * std::sin(phi),
            cos_theta
        };

        std::array<float, 3> pos{
            position.at(0) + r * direction[0],
            position.at(1) + r * direction[1],
            position.at(2) + r * direction[2]
        };

        sampling_position_list[i] = std::make_tuple(r, pos);
    }

    return sampling_position_list;
}
