#include "CylinderSampler.hpp"
#include "Constants.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

#include <random>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace {
    using Vec3 = std::array<float,3>;

    inline Vec3 add(const Vec3& a, const Vec3& b) {
        return { a[0]+b[0], a[1]+b[1], a[2]+b[2] };
    }
    inline Vec3 mul(const Vec3& a, float s) {
        return { a[0]*s, a[1]*s, a[2]*s };
    }
    inline float dot(const Vec3& a, const Vec3& b) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    inline Vec3 cross(const Vec3& a, const Vec3& b) {
        return {
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]
        };
    }
    inline float norm(const Vec3& a) {
        return std::sqrt(dot(a,a));
    }
    inline Vec3 normalize(const Vec3& a) {
        float n = norm(a);
        if (n == 0.0f) return {0.f,0.f,0.f};
        return { a[0]/n, a[1]/n, a[2]/n };
    }

    // 從任意軸向量建立一組正交基 (u, v, w)，其中 w = normalized axis
    // u, v 橫跨圓柱截面平面（皆為單位向量）
    inline void orthonormal_basis_from_axis(const Vec3& axis_raw, Vec3& u, Vec3& v, Vec3& w)
    {
        w = normalize(axis_raw);
        if (w == Vec3{0.f,0.f,0.f})
            throw std::invalid_argument("CylinderSampler: axis_vector cannot be zero.");

        // 避免與 w 幾乎平行，選一個參考向量 t
        Vec3 t = (std::abs(w[2]) < 0.9f) ? Vec3{0.f,0.f,1.f} : Vec3{1.f,0.f,0.f};

        u = cross(w, t);
        float nu = norm(u);
        if (nu == 0.0f) {
            // 萬一 t 恰好平行，換另一個
            t = Vec3{0.f,1.f,0.f};
            u = cross(w, t);
            nu = norm(u);
            if (nu == 0.0f)
                throw std::runtime_error("CylinderSampler: failed to build orthonormal basis.");
        }
        u = mul(u, 1.0f / nu);    // normalize u
        v = cross(w, u);          // v 自然是單位向量且與 u,w 正交
    }
}

CylinderSampler::CylinderSampler(void)
: m_sampling_size{ 10 },
  m_radius_min{ 0.0 },
  m_radius_max{ 1.0 },
  m_height_min{ -0.5 },
  m_height_max{  0.5 }
{
}

void CylinderSampler::Print(void) const
{
    Logger::Log(LogLevel::Info,
        "CylinderSampler Configuration:\n"
        " - Sampling size: " + std::to_string(m_sampling_size) + "\n"
        " - Radius range: ["
        + StringHelper::ToStringWithPrecision<double>(m_radius_min, 3) + ", "
        + StringHelper::ToStringWithPrecision<double>(m_radius_max, 3) + "] Angstrom\n"
        " - Height range (relative): ["
        + StringHelper::ToStringWithPrecision<double>(m_height_min, 3) + ", "
        + StringHelper::ToStringWithPrecision<double>(m_height_max, 3) + "] Angstrom");
}

std::vector<std::tuple<float, std::array<float, 3>>> CylinderSampler::GenerateSamplingPoints(
    const std::array<float, 3> & reference_position,
    const std::array<float, 3> & axis_vector) const
{
    if (m_radius_min < 0.0 || m_radius_max < 0.0)
        throw std::invalid_argument("CylinderSampler: radius range cannot be negative.");
    if (m_radius_min > m_radius_max)
        throw std::invalid_argument("CylinderSampler: radius minimum greater than maximum.");
    if (m_height_min > m_height_max)
        throw std::invalid_argument("CylinderSampler: height minimum greater than maximum.");

    std::array<float,3> u, v, w;
    orthonormal_basis_from_axis(axis_vector, u, v, w);

    std::vector<std::tuple<float, std::array<float, 3>>> out;
    out.resize(m_sampling_size);

    static thread_local std::mt19937 engine{ std::random_device{}() };
    std::uniform_real_distribution<float> dist_r   (static_cast<float>(m_radius_min),
                                                    static_cast<float>(m_radius_max));
    std::uniform_real_distribution<float> dist_phi (0.0f, static_cast<float>(Constants::two_pi));
    std::uniform_real_distribution<float> dist_h   (static_cast<float>(m_height_min),
                                                    static_cast<float>(m_height_max));

    for (unsigned int i = 0; i < m_sampling_size; ++i)
    {
        float r   { dist_r(engine) };
        float phi { dist_phi(engine) };
        float h   { dist_h(engine) };

        std::array<float,3> pos = add(
            add(reference_position, mul(w, h)),
            add(mul(u, r * std::cos(phi)),
                mul(v, r * std::sin(phi)))
        );

        out[i] = std::make_tuple(r, pos);
    }

    return out;
}
