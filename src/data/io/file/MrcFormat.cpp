#include "MrcFormat.hpp"

#include "MapHelper.hpp"

#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <utility>

namespace rhbm_gem {

MrcFormat::MrcFormat()
{
    InitHeader();
}

std::unique_ptr<MapObject> MrcFormat::ReadMap(
    std::istream & stream,
    const std::string & source_name)
{
    InitHeader();
    if (!stream)
    {
        throw std::runtime_error("MrcFormat::Read() failed: invalid input stream.");
    }

    try
    {
        LoadHeader(stream);
        if (m_header.extra_size < 0)
        {
            throw std::runtime_error("MRC extended header size cannot be negative.");
        }
        const std::array<int, 3> array_size{
            m_header.array_size[0],
            m_header.array_size[1],
            m_header.array_size[2]
        };
        const std::array<int, 3> axis_order{
            m_header.axis[0],
            m_header.axis[1],
            m_header.axis[2]
        };
        const auto data_offset{
            static_cast<std::streamoff>(HEAD::SIZE_HEADER)
            + static_cast<std::streamoff>(m_header.extra_size)
        };
        auto data{ data_internal::ReadFloat32VoxelData(
            stream, data_offset, array_size, axis_order, m_header.mode) };
        data_internal::ReorderHeaderAxesToCanonical(
            m_header.array_size,
            m_header.location_index,
            m_header.grid_size,
            m_header.cell_dimension,
            m_header.cell_angle,
            m_header.axis,
            m_header.origin);
        return std::make_unique<MapObject>(
            GetGridSize(), GetGridSpacing(), GetOrigin(), std::move(data));
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "MrcFormat::Read() failed for '" + source_name + "': " + ex.what());
    }
}

void MrcFormat::WriteMap(const MapObject & map_object, std::ostream & stream)
{
    if (!stream)
    {
        throw std::runtime_error("MrcFormat::Write() failed: invalid output stream.");
    }

    InitHeader();
    SetHeader(
        map_object.GetGridSize(),
        map_object.GetGridSpacing(),
        map_object.GetOrigin());
    SaveHeader(stream);
    data_internal::WriteFloat32VoxelData(
        stream,
        static_cast<std::streamoff>(HEAD::SIZE_HEADER),
        map_object.GetMapValueArray(),
        map_object.GetMapValueArraySize(),
        map_object.GetGridSize(),
        m_header.mode);
}

void MrcFormat::InitHeader()
{
    std::memset(&m_header, 0, sizeof(m_header));
    std::fill_n(m_header.array_size, 3, 1);
    m_header.mode = kFloat32Mode;
    std::fill_n(m_header.location_index, 3, 0);
    std::fill_n(m_header.grid_size, 3, 1);
    std::fill_n(m_header.cell_dimension, 3, 1.0f);
    std::fill_n(m_header.cell_angle, 3, 90.0f);
    m_header.axis[0] = 1;
    m_header.axis[1] = 2;
    m_header.axis[2] = 3;
    m_header.map_format_id[0] = 'M';
    m_header.map_format_id[1] = 'A';
    m_header.map_format_id[2] = 'P';
}

void MrcFormat::LoadHeader(std::istream & stream)
{
    stream.seekg(0, std::ios::beg);
    stream.read(reinterpret_cast<char *>(&m_header), sizeof(m_header));
    if (!stream)
    {
        throw std::runtime_error("Failed to read MRC header.");
    }
}

void MrcFormat::SaveHeader(std::ostream & stream)
{
    stream.seekp(0, std::ios::beg);
    stream.write(reinterpret_cast<const char *>(&m_header), sizeof(m_header));
    if (!stream)
    {
        throw std::runtime_error("Failed to write MRC header.");
    }
}

std::array<int, 3> MrcFormat::GetGridSize() const
{
    return {
        m_header.array_size[0],
        m_header.array_size[1],
        m_header.array_size[2]
    };
}

std::array<double, 3> MrcFormat::GetGridSpacing() const
{
    if (m_header.grid_size[0] == 0
        || m_header.grid_size[1] == 0
        || m_header.grid_size[2] == 0)
    {
        throw std::runtime_error("MRC grid size contains a zero dimension.");
    }
    return {
        static_cast<double>(m_header.cell_dimension[0]) / static_cast<double>(m_header.grid_size[0]),
        static_cast<double>(m_header.cell_dimension[1]) / static_cast<double>(m_header.grid_size[1]),
        static_cast<double>(m_header.cell_dimension[2]) / static_cast<double>(m_header.grid_size[2])
    };
}

std::array<double, 3> MrcFormat::GetOrigin() const
{
    return { m_header.origin[0], m_header.origin[1], m_header.origin[2] };
}

void MrcFormat::SetHeader(
    const std::array<int, 3> & grid_size,
    const std::array<double, 3> & grid_spacing,
    const std::array<double, 3> & origin)
{
    numeric_validation::RequireAllPositive(grid_size, "grid_size");
    numeric_validation::RequireAllFinitePositive(grid_spacing, "grid_spacing");

    std::copy(grid_size.begin(), grid_size.end(), m_header.array_size);
    std::copy(grid_size.begin(), grid_size.end(), m_header.grid_size);
    for (size_t i = 0; i < grid_size.size(); ++i)
    {
        m_header.cell_dimension[i] = static_cast<float>(
            grid_spacing[i] * static_cast<double>(grid_size[i]));
        m_header.origin[i] = static_cast<float>(origin[i]);
    }
}

} // namespace rhbm_gem
