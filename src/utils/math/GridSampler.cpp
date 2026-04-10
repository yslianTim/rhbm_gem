#include <rhbm_gem/utils/math/GridSampler.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <stdexcept>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <Eigen/Dense>

using Eigen::Vector3f;

GridSampler::GridSampler() :
    m_sampling_size{ 10 },
    m_window_size{ 5.0 },
    m_reference_u_vector{ 1.0, 0.0, 0.0 }
{
}

void GridSampler::Print() const
{
    std::ostringstream oss;
    oss << "GridSampler Configuration:\n"
        << " - Sampling size: " << m_sampling_size << '\n'
        << " - Window size: "
        << StringHelper::ToStringWithPrecision<double>(m_window_size, 1) << " Angstrom\n";
    Logger::Log(LogLevel::Info, oss.str());
}

SamplingPointList GridSampler::GenerateSamplingPoints(
    const std::array<float, 3> & reference_position,
    const std::array<float, 3> & plane_normal) const
{
    if (m_window_size <= 0.0)
    {
        throw std::invalid_argument("GridSampler: window size should be positive.");
    }
    if (m_sampling_size < 2)
    {
        throw std::invalid_argument("GridSampler: sampling size should be greater than 1.");
    }

    const Eigen::Map<const Vector3f> eigen_reference_position(reference_position.data());
    const Eigen::Map<const Vector3f> eigen_plane_normal(plane_normal.data());
    const Eigen::Map<const Vector3f> eigen_u_vector(m_reference_u_vector.data());

    const auto eps{ Eigen::NumTraits<float>::epsilon() };
    if (eigen_plane_normal.isZero(eps))
    {
        throw std::invalid_argument("GridSampler: plane normal cannot be zero.");
    }
    Vector3f n_unit{ eigen_plane_normal.normalized() };
    Vector3f u_proj{ eigen_u_vector - (eigen_u_vector.dot(n_unit)) * n_unit };
    Vector3f u_unit{ u_proj.normalized() };
    Vector3f v_unit{ n_unit.cross(u_unit) };
    auto half_window_size{ m_window_size / 2.0f };
    auto step_size{ m_window_size / static_cast<float>(m_sampling_size - 1) };

    auto total_grid_size{ m_sampling_size * m_sampling_size };
    SamplingPointList output_list;
    output_list.reserve(total_grid_size);

    for (unsigned int j = 0; j < m_sampling_size; j++)
    {
        float v{ -half_window_size + step_size * static_cast<float>(j) };
        for (unsigned int i = 0; i < m_sampling_size; i++)
        {
            float u{ -half_window_size + step_size * static_cast<float>(i) };
            Vector3f shift{ u * u_unit + v * v_unit };
            Vector3f position{ eigen_reference_position + shift };
            float radius{ shift.norm() };

            output_list.emplace_back(
                radius, std::array<float, 3>{ position.x(), position.y(), position.z() });
        }
    }

    return output_list;
}
