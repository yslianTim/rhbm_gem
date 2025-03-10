#include "SphereSampler.hpp"
#include <Eigen/Dense>
#include <omp.h>

SphereSampler::SphereSampler(void) :
    m_thread_size{ 1 }, m_sampling_size{ 10 }, m_distance_min{ 0.0 }, m_distance_max{ 1.0 },
    m_distance_range{ m_distance_max - m_distance_min }
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

const std::vector<std::tuple<float, std::array<float, 3>>> &  SphereSampler::GenerateSamplingPoints(
    std::array<float, 3> position)
{
    m_sampling_position_list.clear();
    m_sampling_position_list.resize(m_sampling_size);

    Eigen::Array3f atom_location{ position.at(0), position.at(1), position.at(2) };
    Eigen::ArrayXf distance_array{ Eigen::ArrayXf::Zero(m_sampling_size) };
    Eigen::Array3Xf position_array{ Eigen::Array3Xf::Zero(3, m_sampling_size) };
    distance_array.setRandom();
    distance_array = (distance_array + 1.0) * 0.5 * m_distance_range + m_distance_min;
    #pragma omp parallel for num_threads(m_thread_size)
    for (int i = 0; i < m_sampling_size; i++)
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
}