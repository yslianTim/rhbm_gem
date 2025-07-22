#include "MapInterpolationVisitor.hpp"
#include "SphereSampler.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"

#include <algorithm>

MapInterpolationVisitor::MapInterpolationVisitor(SphereSampler * sphere_sampler) :
    m_sphere_sampler{ sphere_sampler }, m_position{ 0.0, 0.0, 0.0 }
{

}

void MapInterpolationVisitor::VisitMapObject(MapObject * data_object)
{
    if (data_object == nullptr) return;
    m_sampling_data_list.clear();
    if (m_sphere_sampler != nullptr)
    {
        const auto & sampling_position_list{ m_sphere_sampler->GenerateSamplingPoints(m_position) };
        m_sampling_data_list.reserve(sampling_position_list.size());
        for (auto & [distance, point] : sampling_position_list)
        {
            auto map_value{ MakeInterpolationInMapObject(data_object, point) };
            m_sampling_data_list.emplace_back(std::make_tuple(distance, map_value));
        }
    }
    else
    {
        m_sampling_data_list.reserve(m_sampling_point_list.size());
        for (auto & [distance, point] : m_sampling_point_list)
        {
            auto map_value{ MakeInterpolationInMapObject(data_object, point) };
            m_sampling_data_list.emplace_back(std::make_tuple(distance, map_value));
        }
    }
}

float MapInterpolationVisitor::MakeInterpolationInMapObject(
    MapObject * data_object, const std::array<float, 3> & position)
{
    auto index{ data_object->GetIndexFromPosition(position) };
    auto origin{ data_object->GetOrigin() };
    auto grid_spacing{ data_object->GetGridSpacing() };

    std::array<float, 3> local;
    for (size_t i = 0; i < 3; i++)
    {
        local.at(i) = (position.at(i) - origin.at(i) - static_cast<float>(index.at(i))*grid_spacing.at(i))/grid_spacing.at(i);
    }

    // Helper function for cubic interpolation
    auto cubic_interpolate = [](float p0, float p1, float p2, float p3, float t)
    {
        float a0 = p1;
        float a1 = 0.5f * (p2 - p0);
        float a2 = 0.5f * (-p3 + 4.0f*p2 - 5.0f*p1 + 2.0f*p0);
        float a3 = 0.5f * (p3 - 3.0f*p2 + 3.0f*p1 - p0);
        return a3 * t * t * t + a2 * t * t + a1 * t + a0;
    };

    // Collect 64 points for interpolation
    std::array<std::array<std::array<float, 4>, 4>, 4> values;
    const auto grid_size{ data_object->GetGridSize() };
    for (int i = -1; i < 3; ++i)
    {
        for (int j = -1; j < 3; ++j)
        {
            for (int k = -1; k < 3; ++k)
            {
                size_t i_next{ static_cast<size_t>(i + 1) };
                size_t j_next{ static_cast<size_t>(j + 1) };
                size_t k_next{ static_cast<size_t>(k + 1) };
                int xi{ std::clamp(index.at(0) + i, 0, grid_size.at(0) - 1) };
                int yj{ std::clamp(index.at(1) + j, 0, grid_size.at(1) - 1) };
                int zk{ std::clamp(index.at(2) + k, 0, grid_size.at(2) - 1) };
                values[i_next][j_next][k_next] = data_object->GetMapValue(xi, yj, zk);
            }
        }
    }
    // Interpolate along x direction
    std::array<std::array<float, 4>, 4> tempY;
    for (size_t j = 0; j < 4; ++j)
    {
        for (size_t k = 0; k < 4; ++k)
        {
            tempY[j][k] = cubic_interpolate(values[0][j][k], values[1][j][k], values[2][j][k], values[3][j][k], local.at(0));
        }
    }

    // Interpolate along y direction
    std::array<float, 4> tempZ;
    for (size_t k = 0; k < 4; ++k)
    {
        tempZ[k] = cubic_interpolate(tempY[0][k], tempY[1][k], tempY[2][k], tempY[3][k], local.at(1));
    }

    // Interpolate along z direction
    return cubic_interpolate(tempZ[0], tempZ[1], tempZ[2], tempZ[3], local.at(2));
}
