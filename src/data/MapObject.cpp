#include "MapObject.hpp"
#include "DataObjectVisitorBase.hpp"
#include "ArrayStats.hpp"

#include <cstring>
#include <cmath>
#include <algorithm>
#include <iomanip>

MapObject::MapObject(
    const std::array<int, 3> & grid_size,
    const std::array<float, 3> & grid_spacing,
    const std::array<float, 3> & origin,
    std::unique_ptr<float[]> map_value_array) :
    m_key_tag{ "" }, m_thread_size{ 4 },
    m_voxel_size{ grid_size.at(0) * grid_size.at(1) * grid_size.at(2) },
    m_map_value_mean{ 0.0f }, m_map_value_min{ 0.0f },
    m_map_value_max{ 0.0f }, m_map_value_sd{ 0.0f },
    m_grid_size{ grid_size }, m_grid_spacing{ grid_spacing }, m_origin{ origin },
    m_map_length{}, m_overflow{}, m_underflow{ origin }, m_upper_bound{}, m_lower_bound{},
    m_map_value_array{ std::move(map_value_array) }
{
    for (int i = 0; i < 3; i++)
    {
        m_map_length.at(i) = static_cast<float>(m_grid_size.at(i) * m_grid_spacing.at(i));
        m_overflow.at(i) = static_cast<float>(m_origin.at(i) + m_map_length.at(i) - m_grid_spacing.at(i));
        m_upper_bound.at(i) = static_cast<float>(m_overflow.at(i) + 0.5 * m_grid_spacing.at(i));
        m_lower_bound.at(i) = static_cast<float>(m_underflow.at(i) - 0.5 * m_grid_spacing.at(i));
    }
    CalculateMapValueMin();
    CalculateMapValueMax();
    CalculateMapValueMean();
    CalculateMapValueSD();
    MapValueArrayNormalization();
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
    CalculateMapValueMin();
    CalculateMapValueMax();
    CalculateMapValueMean();
    CalculateMapValueSD();
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

void MapObject::Accept(DataObjectVisitorBase * visitor)
{
    visitor->VisitMapObject(this);
}

int MapObject::GetGlobalIndex(int index_x, int index_y, int index_z) const
{
    return index_x + m_grid_size.at(0) * (index_y + m_grid_size.at(1) * index_z);
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
    for (int i = 0; i < m_voxel_size; i++)
    {
        m_map_value_array[i] /= m_map_value_sd;
    }
}