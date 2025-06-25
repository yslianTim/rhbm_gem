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

void SphereSampler::Print(void) const
{
    std::cout <<"Sampling size = "<< m_sampling_size << std::endl;
    std::cout <<"Sampling distance from "<< m_distance_min << " to " << m_distance_max <<" A"<< std::endl;
}

std::vector<std::tuple<float, std::array<float, 3>>> SphereSampler::GenerateSamplingPoints(
    const std::array<float, 3> & position)
{
    std::vector<std::tuple<float, std::array<float, 3>>> sampling_position_list;
    sampling_position_list.resize(m_sampling_size);
    auto distance_range{ m_distance_max - m_distance_min };

    Eigen::Array3f atom_location{ position.at(0), position.at(1), position.at(2) };
    Eigen::ArrayXf distance_array{ Eigen::ArrayXf::Zero(m_sampling_size) };
    distance_array.setRandom();
    distance_array = (distance_array + 1.0) * 0.5 * distance_range + m_distance_min;

    #ifdef USE_OPENMP
    #pragma omp parallel for num_threads(m_thread_size)
    #endif
    for (unsigned int i = 0; i < m_sampling_size; i++)
    {
        Eigen::Array3f pos{ Eigen::Quaternionf::UnitRandom() * Eigen::Vector3f::UnitZ() };
        pos *= distance_array(i);
        pos += atom_location;
        sampling_position_list[i] =
            std::make_tuple(distance_array(i), std::array<float, 3>{ pos(0), pos(1), pos(2) });
    }

    return sampling_position_list;
}
