#include "CylinderSampler.hpp"
#include "Constants.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

#include <random>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

using Eigen::Vector3f;

namespace
{
    inline Vector3f ToEigen(const std::array<float, 3> & vec)
    {
        return Vector3f{ vec[0], vec[1], vec[2] };
    }

    inline std::array<float, 3> ToStdArray(const Vector3f & vec)
    {
        return { vec[0], vec[1], vec[2] };
    }

    inline void OrthonormalBasisFromAxis(
        const Vector3f & axis_raw, Vector3f & u, Vector3f & v, Vector3f & w)
    {
        w = axis_raw.normalized();
        if (w == Vector3f::Zero())
        {
            throw std::invalid_argument("CylinderSampler: axis_vector cannot be zero.");
        }

        Vector3f t{ (std::abs(w[2]) < 0.9f) ?
            Vector3f{ 0.0f, 0.0f, 1.0f } : Vector3f{ 1.0f, 0.0f, 0.0f }
        };

        u = w.cross(t);
        float nu{ u.norm() };
        if (nu == 0.0f)
        {
            t = Vector3f{ 0.0f, 1.0f, 0.0f };
            u = w.cross(t);
            nu = u.norm();
            if (nu == 0.0f)
                throw std::runtime_error("CylinderSampler: failed to build orthonormal basis.");
        }
        u /= nu;
        v = w.cross(u);
    }
}

CylinderSampler::CylinderSampler(void) :
    m_sampling_size{ 10 },
    m_distance_min{ 0.0 }, m_distance_max{ 1.0 }, m_height{ 0.1 }
{
}

void CylinderSampler::Print(void) const
{
    Logger::Log(LogLevel::Info,
        "CylinderSampler Configuration:\n"
        " - Sampling size: " + std::to_string(m_sampling_size) + "\n"
        " - Distance range: ["
        + StringHelper::ToStringWithPrecision<double>(m_distance_min, 3) + ", "
        + StringHelper::ToStringWithPrecision<double>(m_distance_max, 3) + "] Angstrom\n"
        " - Height length = "
        + StringHelper::ToStringWithPrecision<double>(m_height, 3) + " Angstrom");
}

std::vector<std::tuple<float, std::array<float, 3>>> CylinderSampler::GenerateSamplingPoints(
    const std::array<float, 3> & reference_position,
    const std::array<float, 3> & axis_vector) const
{
    if (m_distance_min < 0.0 || m_distance_max < 0.0)
    {
        throw std::invalid_argument("CylinderSampler: radius range cannot be negative.");
    }
    if (m_distance_min > m_distance_max)
    {
        throw std::invalid_argument("CylinderSampler: radius minimum greater than maximum.");
    }
    if (m_height < 0.0)
    {
        throw std::invalid_argument("CylinderSampler: height is negative.");
    }

    Vector3f eigen_reference_position{ ToEigen(reference_position) };
    Vector3f eigen_axis_vector{ ToEigen(axis_vector) };

    Vector3f u, v, w;
    OrthonormalBasisFromAxis(eigen_axis_vector, u, v, w);

    std::vector<std::tuple<float, std::array<float, 3>>> output_list;
    output_list.resize(m_sampling_size);

    static thread_local std::mt19937 engine{ std::random_device{}() };
    std::uniform_real_distribution<float> dist_radius(
        static_cast<float>(m_distance_min), static_cast<float>(m_distance_max));
    std::uniform_real_distribution<float> dist_phi(
        0.0f, static_cast<float>(Constants::two_pi));
    std::uniform_real_distribution<float> dist_height(
        static_cast<float>(-0.5 * m_height), static_cast<float>(0.5 * m_height));

    for (unsigned int i = 0; i < m_sampling_size; i++)
    {
        float radius{ dist_radius(engine) };
        float phi{ dist_phi(engine) };
        float height{ dist_height(engine) };

        Vector3f eigen_position{ (
            (eigen_reference_position + (w * height)) +
            ((u * radius * std::cos(phi)) + (v * radius * std::sin(phi))))
        };

        output_list[i] = std::make_tuple(radius, ToStdArray(eigen_position));
    }

    return output_list;
}
