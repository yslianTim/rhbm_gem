#include "MrcFormat.hpp"
#include "Logger.hpp"
#include "MapObject.hpp"

#include <sstream>
#include <fstream>
#include <cstring>
#include <stdexcept>
#include <algorithm>

namespace rhbm_gem {

MrcFormat::MrcFormat()
{
    InitHeader();
}

void MrcFormat::Read(std::istream & stream, const std::string & source_name)
{
    m_data_array.reset();
    InitHeader();
    if (!stream)
    {
        throw std::runtime_error("MrcFormat::Read() failed: invalid input stream.");
    }

    try
    {
        LoadHeader(stream);
        PrintHeader();
        LoadDataArray(stream);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error("MrcFormat::Read() failed for '" + source_name
                                 + "': " + ex.what());
    }
}

void MrcFormat::Write(const MapObject & map_object, std::ostream & stream)
{
    if (!stream)
    {
        throw std::runtime_error("MrcFormat::Write() failed: invalid output stream.");
    }
    const float * data{ map_object.GetMapValueArray() };
    if (data == nullptr)
    {
        throw std::runtime_error("MrcFormat::Write() failed: map value array is null.");
    }

    InitHeader();
    SetHeader(map_object.GetGridSize(), map_object.GetGridSpacing(), map_object.GetOrigin());
    SaveHeader(stream);
    PrintHeader();
    SaveDataArray(data, map_object.GetMapValueArraySize(), stream);
}

void MrcFormat::InitHeader()
{
    std::memset(&m_header, 0, sizeof(m_header));
    std::fill_n(m_header.array_size,        3, 1);
    m_header.mode = static_cast<int>(MODE::SIGNED_FLOAT32);
    std::fill_n(m_header.location_index,    3, 0);
    std::fill_n(m_header.grid_size,         3, 1);
    std::fill_n(m_header.cell_dimension,    3, 1.0f);
    std::fill_n(m_header.cell_angle,        3, 90.0f);
    m_header.axis[0] = 1;
    m_header.axis[1] = 2;
    m_header.axis[2] = 3;
    m_header.min_density = 0.0f;
    m_header.max_density = 0.0f;
    m_header.mean_density = 0.0f;
    m_header.space_group = 0;
    m_header.extra_size = 0;
    std::fill_n(m_header.extra, HEAD::SIZE_WORD, '\0');
    m_header.imodStamp = 0;
    m_header.imodFlags = 0;
    m_header.idtype = 0;
    m_header.lens = 0;
    m_header.nd1 = 0;
    m_header.nd2 = 0;
    m_header.vd1 = 0;
    m_header.vd2 = 0;
    std::fill_n(m_header.tiltangles, 6, 0.0f);
    m_header.map_format_id[0] = 'M';
    m_header.map_format_id[1] = 'A';
    m_header.map_format_id[2] = 'P';
    m_header.map_format_id[3] = '\0';
    std::fill_n(m_header.machine_stamp, 4, '\0');
    m_header.rms = 0.0f;
    m_header.label_size = 0;
    std::fill_n(&m_header.label[0][0], HEAD::NUM_LABEL * HEAD::SIZE_LABEL, '\0');
}

void MrcFormat::LoadHeader(std::istream & stream)
{
    stream.seekg(0, std::ios::beg);
    stream.read(reinterpret_cast<char*>(&m_header), sizeof(m_header));
    if (!stream)
    {
        throw std::runtime_error("LoadHeader failed!");
    }
}

void MrcFormat::SaveHeader(std::ostream & stream)
{
    stream.seekp(0, std::ios::beg);
    stream.write(reinterpret_cast<const char*>(&m_header), sizeof(m_header));
    if (!stream)
    {
        throw std::runtime_error("SaveHeader failed!");
    }
}

void MrcFormat::PrintHeader() const
{
    Logger::Log(LogLevel::Debug,
        "MRC Header Information:\n"
        "Array Size: " + std::to_string(m_header.array_size[0]) + " x "
        + std::to_string(m_header.array_size[1]) + " x "
        + std::to_string(m_header.array_size[2]) + "\n"
        "Mode: " + std::to_string(m_header.mode) + "\n"
        "Location Index: " + std::to_string(m_header.location_index[0]) + " x "
        + std::to_string(m_header.location_index[1]) + " x "
        + std::to_string(m_header.location_index[2]) + "\n"
        "Grid Size: " + std::to_string(m_header.grid_size[0]) + " x "
        + std::to_string(m_header.grid_size[1]) + " x "
        + std::to_string(m_header.grid_size[2]) + "\n"
        "Map Length: " + std::to_string(m_header.cell_dimension[0]) + ", "
        + std::to_string(m_header.cell_dimension[1]) + ", "
        + std::to_string(m_header.cell_dimension[2]) + "\n"
        "Cell Angles: " + std::to_string(m_header.cell_angle[0]) + ", "
        + std::to_string(m_header.cell_angle[1]) + ", "
        + std::to_string(m_header.cell_angle[2]) + "\n"
        "Axis: " + std::to_string(m_header.axis[0]) + ", "
        + std::to_string(m_header.axis[1]) + ", "
        + std::to_string(m_header.axis[2]) + "\n"
        "Min Density: " + std::to_string(m_header.min_density) + "\n"
        "Max Density: " + std::to_string(m_header.max_density) + "\n"
        "Mean Density: " + std::to_string(m_header.mean_density) + "\n"
        "Space Group: " + std::to_string(m_header.space_group) + "\n"
        "Extra Size: " + std::to_string(m_header.extra_size) + "\n"
        "Tilt angles: " + std::to_string(m_header.tiltangles[0]) + ", "
        + std::to_string(m_header.tiltangles[1]) + ", "
        + std::to_string(m_header.tiltangles[2]) + ", "
        + std::to_string(m_header.tiltangles[3]) + ", "
        + std::to_string(m_header.tiltangles[4]) + ", "
        + std::to_string(m_header.tiltangles[5]) + "\n"
        "Origin: " + std::to_string(m_header.origin[0]) + ", "
        + std::to_string(m_header.origin[1]) + ", "
        + std::to_string(m_header.origin[2]) + "\n"
        "Map Format ID: " + std::string(m_header.map_format_id) + "\n"
        "Machine Stamp: " + std::string(m_header.machine_stamp) + "\n"
        "RMS: " + std::to_string(m_header.rms) + "\n"
        "Label Size: " + std::to_string(m_header.label_size) + "\n"
    );
}

size_t MrcFormat::GetElementSize() const
{
    switch (static_cast<MODE>(m_header.mode))
    {
        case MODE::SIGNED_INT8:
            return 1;
        case MODE::SIGNED_INT16:
        case MODE::UNSIGNED_INT16:
        case MODE::IEEE754_FLOAT16:
            return 2;
        case MODE::SIGNED_FLOAT32:
            return 4;
        case MODE::COMPLEX_INT16:
            return 4;
        case MODE::COMPLEX_FLOAT32:
            return 8;
        case MODE::TWOPERBYTE_BIT4:
            throw std::runtime_error("MODE TWOPERBYTE_BIT4 is not supported.");
        default:
            throw std::runtime_error("Unknown mode value");
    }
}

void MrcFormat::LoadDataArray(std::istream & stream)
{
    // Ensure we start reading data from the beginning of the file
    stream.seekg(HEAD::SIZE_HEADER, std::ios::beg);
    
    // Skip if extra header is exist（extra_size > 0)
    if (m_header.extra_size > 0)
    {
        stream.seekg(m_header.extra_size, std::ios::cur);
    }
    
    size_t num_voxels{
        static_cast<size_t>(m_header.array_size[0]) *
        static_cast<size_t>(m_header.array_size[1]) *
        static_cast<size_t>(m_header.array_size[2])
    };
    
    size_t element_size{ GetElementSize() };
    size_t total_bytes{ num_voxels * element_size };
    
    // Read raw data into temporary buffer first
    auto raw_data{ std::make_unique<float[]>(num_voxels) };
    switch (static_cast<MODE>(m_header.mode))
    {
        case MODE::SIGNED_FLOAT32:
            stream.read(reinterpret_cast<char*>(raw_data.get()),
                        static_cast<std::streamsize>(total_bytes));
            if (!stream)
            {
                throw std::runtime_error("Failed to read voxel data from file");
            }
            break;
        default:
            throw std::runtime_error("Unsupported MODE in LoadDataArray");
    }

    if (m_header.axis[0] == 1 && m_header.axis[1] == 2 && m_header.axis[2] == 3)
    {
        // Data is already in X->Y->Z order, no reordering needed
        m_data_array = std::move(raw_data);
        return;
    }

    // Build mapping from X/Y/Z axis to column/row/section positions
    int axis_to_index[3];
    for (int i = 0; i < 3; i++)
    {
        // axis values are 1-based, convert to 0-based indices
        axis_to_index[m_header.axis[i] - 1] = i;
    }

    // Determine dimension sizes in canonical X,Y,Z order
    size_t dims[3]{
        static_cast<size_t>(m_header.array_size[axis_to_index[0]]),
        static_cast<size_t>(m_header.array_size[axis_to_index[1]]),
        static_cast<size_t>(m_header.array_size[axis_to_index[2]])
    };

    // Compute strides for each axis in the source buffer
    size_t src_stride[3];
    size_t stride_acc{ 1 };
    for (int i = 0; i < 3; ++i)
    {
        src_stride[m_header.axis[i] - 1] = stride_acc;
        stride_acc *= static_cast<size_t>(m_header.array_size[i]);
    }

    // Allocate destination array and reorder data into X->Y->Z order
    auto reordered_array{ std::make_unique<float[]>(num_voxels) };
    for (size_t z = 0; z < dims[2]; z++)
    {
        size_t src_off_z{ z * src_stride[2] };
        size_t dst_off_z{ z * dims[0] * dims[1] };
        for (size_t y = 0; y < dims[1]; y++)
        {
            size_t src_off_y{ src_off_z + y * src_stride[1] };
            size_t dst_off_y{ dst_off_z + y * dims[0] };
            for (size_t x = 0; x < dims[0]; x++)
            {
                reordered_array[dst_off_y + x] = raw_data[src_off_y + x * src_stride[0]];
            }
        }
    }

    // Update header to reflect canonical axis order and dimensions
    ReorderedAxisRelatedParameters();
    
    m_data_array = std::move(reordered_array);
}

void MrcFormat::SaveDataArray(const float * data, size_t size, std::ostream & stream)
{
    size_t expected_voxels{ static_cast<size_t>(m_header.array_size[0]) *
                            static_cast<size_t>(m_header.array_size[1]) *
                            static_cast<size_t>(m_header.array_size[2]) };
    if (size != expected_voxels)
    {
        throw std::runtime_error("SaveDataArray: voxel count does not match header");
    }

    std::streamoff data_offset{ HEAD::SIZE_HEADER + static_cast<std::streamoff>(m_header.extra_size) };
    stream.seekp(data_offset, std::ios::beg);
    if (!stream)
    {
        throw std::runtime_error("Failed to seek to data offset");
    }

    const size_t element_size{ GetElementSize() };
    const size_t total_bytes{ size * element_size };

    stream.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(total_bytes));
    if (!stream)
    {
        throw std::runtime_error("Failed to write data array to file");
    }
}

std::unique_ptr<float[]> MrcFormat::GetDataArray()
{
    return std::move(m_data_array);
}

std::array<int, 3> MrcFormat::GetGridSize()
{
    std::array<int, 3> grid_size{
        m_header.array_size[0],
        m_header.array_size[1],
        m_header.array_size[2]
    };
    return grid_size;
}

std::array<float, 3> MrcFormat::GetGridSpacing()
{
    if (m_header.grid_size[0] == 0 || m_header.grid_size[1] == 0 || m_header.grid_size[2] == 0)
    {
        throw std::runtime_error("GetGridSpacing: grid_size has zero dimension");
    }
    std::array<float, 3> grid_spacing{
        m_header.cell_dimension[0] / static_cast<float>(m_header.grid_size[0]),
        m_header.cell_dimension[1] / static_cast<float>(m_header.grid_size[1]),
        m_header.cell_dimension[2] / static_cast<float>(m_header.grid_size[2])
    };
    return grid_spacing;
}

std::array<float, 3> MrcFormat::GetOrigin()
{
    std::array<float, 3> origin{
        m_header.origin[0],
        m_header.origin[1],
        m_header.origin[2]
    };
    return origin;
}

void MrcFormat::SetHeader(const std::array<int, 3> & grid_size,
                          const std::array<float, 3> & grid_spacing,
                          const std::array<float, 3> & origin)
{
    if (grid_size[0] <= 0 || grid_size[1] <= 0 || grid_size[2] <= 0)
    {
        throw std::runtime_error("SetHeader: grid_size must be positive");
    }
    if (grid_spacing[0] <= 0.0f || grid_spacing[1] <= 0.0f || grid_spacing[2] <= 0.0f)
    {
        throw std::runtime_error("SetHeader: grid_spacing must be positive");
    }

    std::memcpy(m_header.array_size, grid_size.data(), sizeof(m_header.array_size));
    std::memcpy(m_header.grid_size, grid_size.data(), sizeof(m_header.grid_size));

    m_header.cell_dimension[0] = grid_spacing[0] * static_cast<float>(m_header.grid_size[0]);
    m_header.cell_dimension[1] = grid_spacing[1] * static_cast<float>(m_header.grid_size[1]);
    m_header.cell_dimension[2] = grid_spacing[2] * static_cast<float>(m_header.grid_size[2]);

    std::memcpy(m_header.origin, origin.data(), sizeof(m_header.origin));
}

void MrcFormat::ReorderedAxisRelatedParameters()
{
    if (m_header.axis[0] == 1 && m_header.axis[1] == 2 && m_header.axis[2] == 3)
    {
        Logger::Log(LogLevel::Debug,
            "MrcFormat::ReorderedAxisRelatedParameters : "
            "Axis already in X->Y->Z order, no need to reorder."
        );
        return;
    }

    int axis_to_index[3];
    for (int i = 0; i < 3; i++)
    {
        // axis values are 1-based, convert to 0-based indices
        axis_to_index[m_header.axis[i] - 1] = i;
    }

    int array_size[3]{
        m_header.array_size[0], m_header.array_size[1], m_header.array_size[2]
    };
    int location_index[3]{
        m_header.location_index[0], m_header.location_index[1], m_header.location_index[2]
    };
    int grid_size[3]{
        m_header.grid_size[0], m_header.grid_size[1], m_header.grid_size[2]
    };
    float cell_dimension[3]{
        m_header.cell_dimension[0], m_header.cell_dimension[1], m_header.cell_dimension[2]
    };
    float cell_angle[3]{
        m_header.cell_angle[0], m_header.cell_angle[1], m_header.cell_angle[2]
    };
    float origin[3]{
        m_header.origin[0], m_header.origin[1], m_header.origin[2]
    };

    // Reorder parameters based on current axis mapping
    for (int i = 0; i < 3; i++)
    {
        m_header.array_size[i] = array_size[axis_to_index[i]];
        m_header.location_index[i] = location_index[axis_to_index[i]];
        m_header.grid_size[i] = grid_size[axis_to_index[i]];
        m_header.cell_dimension[i] = cell_dimension[axis_to_index[i]];
        m_header.cell_angle[i] = cell_angle[axis_to_index[i]];
        m_header.origin[i] = origin[axis_to_index[i]];
    }

    // Update axis to canonical X->Y->Z order
    m_header.axis[0] = 1;
    m_header.axis[1] = 2;
    m_header.axis[2] = 3;
}

} // namespace rhbm_gem
