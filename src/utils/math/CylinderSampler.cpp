#include <rhbm_gem/utils/math/CylinderSampler.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <random>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <Eigen/Dense>

using Eigen::Vector3f;

CylinderSampler::CylinderSampler() :
    m_sampling_size{ 10 },
    m_distance_min{ 0.0 }, m_distance_max{ 1.0 }, m_height{ 0.1 }
{
}

void CylinderSampler::Print() const
{
    std::ostringstream oss;
    oss << "CylinderSampler Configuration:\n"
        << " - Sampling size: " << m_sampling_size << '\n'
        << " - Distance range: ["
        << StringHelper::ToStringWithPrecision<double>(m_distance_min, 1) << ", "
        << StringHelper::ToStringWithPrecision<double>(m_distance_max, 1) << "] Angstrom\n"
        << " - Height length = "
        << StringHelper::ToStringWithPrecision<double>(m_height, 1) << " Angstrom";
    Logger::Log(LogLevel::Info, oss.str());
}

void CylinderSampler::SetDistanceRange(double min_value, double max_value)
{
    if (min_value < 0.0 || max_value < 0.0)
    {
        throw std::invalid_argument("CylinderSampler: radius range cannot be negative.");
    }
    if (min_value > max_value)
    {
        throw std::invalid_argument("CylinderSampler: radius minimum greater than maximum.");
    }

    m_distance_min = min_value;
    m_distance_max = max_value;
}

SamplingPointList CylinderSampler::GenerateSamplingPoints(
    const std::array<float, 3> & reference_position,
    const std::array<float, 3> & axis_vector) const
{
    if (m_height < 0.0)
    {
        throw std::invalid_argument("CylinderSampler: height is negative.");
    }

    const Eigen::Map<const Vector3f> eigen_reference_position(reference_position.data());
    const Eigen::Map<const Vector3f> eigen_axis_vector(axis_vector.data());

    Vector3f u, v, w;
    if (eigen_axis_vector.isZero(Eigen::NumTraits<float>::epsilon()))
    {
        throw std::invalid_argument("CylinderSampler: axis_vector cannot be zero.");
    }
    w = eigen_axis_vector.normalized();
    u = w.unitOrthogonal();
    v = w.cross(u);

    SamplingPointList output_list;
    output_list.reserve(m_sampling_size);

    static thread_local std::mt19937 engine{ std::random_device{}() };
    const float min_radius{ static_cast<float>(m_distance_min) };
    const float max_radius{ static_cast<float>(m_distance_max) };
    std::uniform_real_distribution<float> dist_radius_squared(
        min_radius * min_radius, max_radius * max_radius);
    std::uniform_real_distribution<float> dist_phi(
        0.0f, static_cast<float>(Constants::two_pi));
    std::uniform_real_distribution<float> dist_height(
        static_cast<float>(-0.5 * m_height), static_cast<float>(0.5 * m_height));

    for (unsigned int i = 0; i < m_sampling_size; i++)
    {
        float radius_squared{ dist_radius_squared(engine) };
        float radius{ std::sqrt(radius_squared) };
        float phi{ dist_phi(engine) };
        float height{ dist_height(engine) };

        Vector3f position{ (
            (eigen_reference_position + (w * height)) +
            ((u * radius * std::cos(phi)) + (v * radius * std::sin(phi))))
        };

        output_list.emplace_back(
            radius, std::array<float, 3>{ position.x(), position.y(), position.z() });
    }

    return output_list;
}
