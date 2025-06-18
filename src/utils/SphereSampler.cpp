#include "SphereSampler.hpp"
#include "Constants.hpp"
#include <Eigen/Dense>

#include <iostream>
#include <random>

#ifdef USE_OPENMP
#include <omp.h>
#endif

SphereSampler::SphereSampler(void) :
    m_thread_size{ 1 }, m_sampling_size{ 10 }, m_distance_min{ 0.0 }, m_distance_max{ 1.0 }
{

}

SphereSampler::~SphereSampler()
{

}

void SphereSampler::Print(void) const
{
    std::cout <<"Sampling size = "<< m_sampling_size << std::endl;
    std::cout <<"Sampling distance from "<< m_distance_min << " to " << m_distance_max <<" A"<< std::endl;
}
/*
const std::vector<std::tuple<float, std::array<float, 3>>> &  SphereSampler::GenerateSamplingPoints(
    std::array<float, 3> position)
{
    m_sampling_position_list.clear();
    m_sampling_position_list.resize(m_sampling_size);
    auto distance_range{ m_distance_max - m_distance_min };

    Eigen::Array3f atom_location{ position.at(0), position.at(1), position.at(2) };
    Eigen::ArrayXf distance_array{ Eigen::ArrayXf::Zero(m_sampling_size) };
    Eigen::Array3Xf position_array{ Eigen::Array3Xf::Zero(3, m_sampling_size) };
    distance_array.setRandom();
    distance_array = (distance_array + 1.0) * 0.5 * distance_range + m_distance_min;

    #ifdef USE_OPENMP
    #pragma omp parallel for num_threads(m_thread_size)
    #endif
    for (unsigned int i = 0; i < m_sampling_size; i++)
    {
        position_array.col(i) = Eigen::Quaternionf::UnitRandom() * Eigen::Vector3f::UnitZ();
        position_array.col(i) *= distance_array(i);
        position_array.col(i) += atom_location;
        m_sampling_position_list.at(i) =
            std::make_tuple(
                distance_array(i),
                std::array<float, 3>{
                    position_array.col(i)(0),
                    position_array.col(i)(1),
                    position_array.col(i)(2)
                });
    }

    return m_sampling_position_list;
}*/

std::vector<std::tuple<float, std::array<float, 3>>>
SphereSampler::GenerateSamplingPoints(std::array<float,3> position)
{
    const float range{ static_cast<float>(m_distance_max - m_distance_min) };
    const Eigen::Vector3f center{ position[0], position[1], position[2] };

    std::vector<std::tuple<float, std::array<float, 3>>> samples;
    samples.resize(m_sampling_size);

    #pragma omp parallel
    {
        thread_local std::mt19937 rng{std::random_device{}()};
        std::uniform_real_distribution<float> U(-1.f, 1.f);
        std::uniform_real_distribution<float> U01(0.f, 1.f);

        // 把平行區工作分配給各執行緒
        #pragma omp for schedule(static)
        for (unsigned i = 0; i < m_sampling_size; ++i)
        {
            float d = U01(rng) * range + static_cast<float>(m_distance_min);

            // Marsaglia 球面均勻取樣
            float  z = U(rng);
            float  t = U(rng);
            float  r = std::sqrt(1.f - z*z);
            float  x = r * std::cos(t * static_cast<float>(Constants::pi));
            float  y = r * std::sin(t * static_cast<float>(Constants::pi));

            std::array<float,3> p{
                center.x() + d * x,
                center.y() + d * y,
                center.z() + d * z };

            samples[i] = std::make_tuple(d, p);
        }
    }
    return samples;
}
