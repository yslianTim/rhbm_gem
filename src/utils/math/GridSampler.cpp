#include <rhbm_gem/utils/math/GridSampler.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <stdexcept>
#include <sstream>
#include <Eigen/Dense>

using Eigen::Vector3d;

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
        << rhbm_gem::string_helper::ToStringWithPrecision<double>(m_window_size, 1) << " Angstrom\n";
    Logger::Log(LogLevel::Info, oss.str());
}

SamplingPointList GridSampler::GenerateSamplingPoints(
    const std::array<double, 3> & reference_position,
    const std::array<double, 3> & plane_normal) const
{
    rhbm_gem::numeric_validation::RequireFinitePositive(
        m_window_size,
        "GridSampler window size");
    rhbm_gem::numeric_validation::RequireAtLeast(
        m_sampling_size,
        2u,
        "GridSampler sampling size");

    const Eigen::Map<const Vector3d> eigen_reference_position(reference_position.data());
    const Eigen::Map<const Vector3d> eigen_plane_normal(plane_normal.data());
    const Eigen::Map<const Vector3d> eigen_u_vector(m_reference_u_vector.data());

    const auto eps{ Eigen::NumTraits<double>::epsilon() };
    if (eigen_plane_normal.isZero(eps))
    {
        throw std::invalid_argument("GridSampler: plane normal cannot be zero.");
    }
    Vector3d n_unit{ eigen_plane_normal.normalized() };
    Vector3d u_proj{ eigen_u_vector - (eigen_u_vector.dot(n_unit)) * n_unit };
    Vector3d u_unit{ u_proj.normalized() };
    Vector3d v_unit{ n_unit.cross(u_unit) };
    auto half_window_size{ m_window_size / 2.0 };
    auto step_size{ m_window_size / static_cast<double>(m_sampling_size - 1) };

    auto total_grid_size{ m_sampling_size * m_sampling_size };
    SamplingPointList output_list;
    output_list.reserve(total_grid_size);

    for (unsigned int j = 0; j < m_sampling_size; j++)
    {
        double v{ -half_window_size + step_size * static_cast<double>(j) };
        for (unsigned int i = 0; i < m_sampling_size; i++)
        {
            double u{ -half_window_size + step_size * static_cast<double>(i) };
            Vector3d shift{ u * u_unit + v * v_unit };
            Vector3d position{ eigen_reference_position + shift };
            double radius{ shift.norm() };

            output_list.emplace_back(SamplingPoint{
                radius,
                std::array<double, 3>{ position.x(), position.y(), position.z() }
            });
        }
    }

    return output_list;
}
