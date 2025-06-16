#include "MapObject.hpp"
#include "DataObjectVisitorBase.hpp"
#include "ArrayStats.hpp"

#include <iostream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <iomanip>

MapObject::MapObject(
    const std::array<int, 3> & grid_size,
    const std::array<float, 3> & grid_spacing,
    const std::array<float, 3> & origin) :
    m_key_tag{ "" }, m_thread_size{ 4 },
    m_voxel_size{ static_cast<size_t>(grid_size.at(0) * grid_size.at(1) * grid_size.at(2)) },
    m_map_value_mean{ 0.0f }, m_map_value_min{ 0.0f },
    m_map_value_max{ 0.0f }, m_map_value_sd{ 0.0f },
    m_grid_size{ grid_size }, m_grid_spacing{ grid_spacing }, m_origin{ origin },
    m_map_length{}, m_overflow{}, m_underflow{ origin }, m_upper_bound{}, m_lower_bound{},
    m_map_value_array{ std::make_unique<float[]>(m_voxel_size) }
{
    for (size_t i = 0; i < 3; i++)
    {
        m_map_length.at(i) = static_cast<float>(m_grid_size.at(i)) * m_grid_spacing.at(i);
        m_overflow.at(i) = static_cast<float>(m_origin.at(i) + m_map_length.at(i) - m_grid_spacing.at(i));
        m_upper_bound.at(i) = static_cast<float>(m_overflow.at(i) + 0.5 * m_grid_spacing.at(i));
        m_lower_bound.at(i) = static_cast<float>(m_underflow.at(i) - 0.5 * m_grid_spacing.at(i));
    }
}

MapObject::MapObject(
    const std::array<int, 3> & grid_size,
    const std::array<float, 3> & grid_spacing,
    const std::array<float, 3> & origin,
    std::unique_ptr<float[]> map_value_array) :
    m_key_tag{ "" }, m_thread_size{ 4 },
    m_voxel_size{ static_cast<size_t>(grid_size.at(0) * grid_size.at(1) * grid_size.at(2)) },
    m_map_value_mean{ 0.0f }, m_map_value_min{ 0.0f },
    m_map_value_max{ 0.0f }, m_map_value_sd{ 0.0f },
    m_grid_size{ grid_size }, m_grid_spacing{ grid_spacing }, m_origin{ origin },
    m_map_length{}, m_overflow{}, m_underflow{ origin }, m_upper_bound{}, m_lower_bound{},
    m_map_value_array{ std::move(map_value_array) }
{
    for (size_t i = 0; i < 3; i++)
    {
        m_map_length.at(i) = static_cast<float>(m_grid_size.at(i)) * m_grid_spacing.at(i);
        m_overflow.at(i) = static_cast<float>(m_origin.at(i) + m_map_length.at(i) - m_grid_spacing.at(i));
        m_upper_bound.at(i) = static_cast<float>(m_overflow.at(i) + 0.5 * m_grid_spacing.at(i));
        m_lower_bound.at(i) = static_cast<float>(m_underflow.at(i) - 0.5 * m_grid_spacing.at(i));
    }
    Update();
}

MapObject::~MapObject()
{

}

MapObject::MapObject(const MapObject & other) :
    m_key_tag{ other.m_key_tag },
    m_thread_size{ other.m_thread_size }, m_voxel_size{ other.m_voxel_size },
    m_map_value_mean{ other.m_map_value_mean }, m_map_value_min{ other.m_map_value_min },
    m_map_value_max{ other.m_map_value_max }, m_map_value_sd{ other.m_map_value_sd },
    m_grid_size{ other.m_grid_size }, m_grid_spacing{ other.m_grid_spacing }, m_origin{ other.m_origin },
    m_map_length{ other.m_map_length }, m_overflow{ other.m_overflow }, m_underflow{ other.m_underflow },
    m_upper_bound{ other.m_upper_bound }, m_lower_bound{ other.m_lower_bound },
    m_map_value_array{ std::make_unique<float[]>(other.m_voxel_size) }
{
    std::memcpy(m_map_value_array.get(), other.m_map_value_array.get(), m_voxel_size * sizeof(float));
    Update();
}

std::unique_ptr<DataObjectBase> MapObject::Clone() const
{
    return std::make_unique<MapObject>(*this);
}

void MapObject::Display(void) const
{
    std::cout <<" o=====================================================o\n";
    std::cout <<" |  Map Object  |   X-axis   |   Y-axis   |   Z-axis   |\n";
    std::cout <<" o=====================================================o\n";
    std::cout <<" | Grid size    | ";
    std::cout << std::setw(10) << m_grid_size.at(0) <<" | "
              << std::setw(10) << m_grid_size.at(1) <<" | "
              << std::setw(10) << m_grid_size.at(2) <<" |\n";
    std::cout <<" | Grid Spacing | ";
    std::cout << std::setw(10) << m_grid_spacing.at(0) <<" | "
              << std::setw(10) << m_grid_spacing.at(1) <<" | "
              << std::setw(10) << m_grid_spacing.at(2) <<" |\n";
    std::cout <<" | Origin (A)   | ";
    std::cout << std::setw(10) << m_origin.at(0) <<" | "
              << std::setw(10) << m_origin.at(1) <<" | "
              << std::setw(10) << m_origin.at(2) <<" |\n";
    std::cout <<" | Map Length(A)| ";
    std::cout << std::setw(10) << m_map_length.at(0) <<" | "
              << std::setw(10) << m_map_length.at(1) <<" | "
              << std::setw(10) << m_map_length.at(2) <<" |\n";
    std::cout <<" |-----------------------------------------------------|\n";
    std::cout <<" | Map value min  | "<< std::setw(34) << m_map_value_min <<" |\n";
    std::cout <<" | Map value max  | "<< std::setw(34) << m_map_value_max <<" |\n";
    std::cout <<" | Map value mean | "<< std::setw(34) << m_map_value_mean <<" |\n";
    std::cout <<" | Map value s.d. | "<< std::setw(34) << m_map_value_sd <<" |\n";
    std::cout <<" o=====================================================o\n";
}

void MapObject::Update(void)
{
    CalculateMapValueMin();
    CalculateMapValueMax();
    CalculateMapValueMean();
    CalculateMapValueSD();
}

void MapObject::Accept(DataObjectVisitorBase * visitor)
{
    visitor->VisitMapObject(this);
}

void MapObject::SetMapValueArray(std::unique_ptr<float[]> map_value_array)
{
    if (m_map_value_array != nullptr && m_map_value_mean != 0.0f)
    {
        std::cout <<"[Warning] MapObject::SetMapValueArray - MapObject already has a map value array."
                  <<" The map value array will be replaced by new inserted data."<< std::endl;
    }
    m_map_value_array = std::move(map_value_array);
    Update();
}

size_t MapObject::GetGlobalIndex(int index_x, int index_y, int index_z) const
{
    if (index_x < 0 || index_x >= m_grid_size[0] ||
        index_y < 0 || index_y >= m_grid_size[1] ||
        index_z < 0 || index_z >= m_grid_size[2])
    {
        throw std::out_of_range(
            "GetGlobalIndex: ("
            + std::to_string(index_x)+","
            + std::to_string(index_y)+","
            + std::to_string(index_z)+") out of grid range!");
    }
    return static_cast<size_t>(index_x + m_grid_size.at(0) * (index_y + m_grid_size.at(1) * index_z));
}

std::array<int, 3> MapObject::GetGridIndex(size_t global_index) const
{
    if (global_index >= m_voxel_size)
    {
        throw std::out_of_range(
            "GetGridIndex: global_index = " + std::to_string(global_index) 
            +" is out of range [0," + std::to_string(m_voxel_size - 1) + "]");
    }
    size_t plane_size{ static_cast<size_t>(m_grid_size[0]) * static_cast<size_t>(m_grid_size[1]) };
    auto z{ static_cast<int>(global_index / plane_size) };
    auto rem{ static_cast<int>(global_index % plane_size) };
    auto y{ rem / m_grid_size[0] };
    auto x{ rem % m_grid_size[0] };
    return std::array{ x, y, z };
}

std::array<float, 3> MapObject::GetGridPosition(size_t global_index) const
{
    auto grid_index{ GetGridIndex(global_index) };
    auto x_pos{ static_cast<float>(grid_index[0]) * m_grid_spacing[0] + m_origin[0] };
    auto y_pos{ static_cast<float>(grid_index[1]) * m_grid_spacing[1] + m_origin[1] };
    auto z_pos{ static_cast<float>(grid_index[2]) * m_grid_spacing[2] + m_origin[2] };
    return std::array{ x_pos, y_pos, z_pos };
}

std::array<int, 3> MapObject::GetIndexFromPosition(const std::array<float, 3> & position) const
{
    try
    {
        CheckPosition(position);
    }
    catch (const std::out_of_range & exception)
    {
        std::cerr << exception.what() << std::endl;
        return std::array<int, 3>{-1, -1, -1};
    }
    auto x{ (position.at(0) - m_origin.at(0)) / m_grid_spacing.at(0) };
    auto y{ (position.at(1) - m_origin.at(1)) / m_grid_spacing.at(1) };
    auto z{ (position.at(2) - m_origin.at(2)) / m_grid_spacing.at(2) };
    return std::array<int, 3>{
        static_cast<int>(std::floor(x)),
        static_cast<int>(std::floor(y)),
        static_cast<int>(std::floor(z))};
}

float MapObject::GetMapValue(size_t global_index) const
{
    return m_map_value_array[global_index];
}

float MapObject::GetMapValue(int index_x, int index_y, int index_z) const
{
    try
    {
        CheckIndex(index_x, index_y, index_z);
    }
    catch(const std::out_of_range & exception)
    {
        std::cerr << exception.what() << '\n';
        return 0.0f;
    }
    return m_map_value_array[GetGlobalIndex(index_x, index_y, index_z)];
}

void MapObject::CheckIndex(int index_x, int index_y, int index_z) const
{
    if (index_x < 0 || index_y < 0 || index_z < 0)
    {
        throw std::out_of_range("Grid index is not positive value.");
    }
    if (index_x >= m_grid_size.at(0) || index_y >= m_grid_size.at(1) || index_z >= m_grid_size.at(2))
    {
        throw std::out_of_range("Grid index is out of array size.");
    }
}

void MapObject::CheckPosition(const std::array<float, 3> & position) const
{
    if (position.at(0) > m_upper_bound.at(0) ||
        position.at(1) > m_upper_bound.at(1) ||
        position.at(2) > m_upper_bound.at(2) ||
        position.at(0) < m_lower_bound.at(0) ||
        position.at(1) < m_lower_bound.at(1) ||
        position.at(2) < m_lower_bound.at(2))
    {
        throw std::out_of_range("Position is out of map boundary.");
    }
}

void MapObject::CalculateMapValueMean(void)
{
    m_map_value_mean = ArrayStats<float>::ComputeMean(
        m_map_value_array.get(), m_voxel_size, m_thread_size);
}

void MapObject::CalculateMapValueMin(void)
{
    m_map_value_min = ArrayStats<float>::ComputeMin(
        m_map_value_array.get(), m_voxel_size, m_thread_size);
}

void MapObject::CalculateMapValueMax(void)
{
    m_map_value_max = ArrayStats<float>::ComputeMax(
        m_map_value_array.get(), m_voxel_size, m_thread_size);
}

void MapObject::CalculateMapValueSD(void)
{
    m_map_value_sd = ArrayStats<float>::ComputeStandardDeviation(
        m_map_value_array.get(), m_voxel_size, m_map_value_mean, m_thread_size);
}

void MapObject::MapValueArrayNormalization(void)
{
    for (size_t i = 0; i < static_cast<size_t>(m_voxel_size); i++)
    {
        m_map_value_array[i] /= m_map_value_sd;
    }
    Update();
}