#include "MrcFormat.hpp"
#include "Logger.hpp"

#include <sstream>
#include <fstream>
#include <cstring>
#include <stdexcept>
#include <algorithm>

MrcFormat::MrcFormat(void)
{
    InitHeader();
}

void MrcFormat::InitHeader(void)
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
    m_header.map[0] = 'M';
    m_header.map[1] = 'A';
    m_header.map[2] = 'P';
    m_header.map[3] = '\0';
    std::fill_n(m_header.stamp, 4, '\0');
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

void MrcFormat::PrintHeader(void) const
{
    std::ostringstream oss;
    oss << "MRC Header Information:\n";
    oss << "Grid Size: " << m_header.grid_size[0] << " x "
              << m_header.grid_size[1] << " x " << m_header.grid_size[2] <<"\n";
    oss << "Cell Dimensions: " << m_header.cell_dimension[0] << ", "
              << m_header.cell_dimension[1] << ", " << m_header.cell_dimension[2] <<"\n";
    Logger::Log(LogLevel::Info, oss.str());
}

size_t MrcFormat::GetElementSize(void) const
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
    
    size_t num_voxels{ static_cast<size_t>(m_header.array_size[0]) *
                       static_cast<size_t>(m_header.array_size[1]) *
                       static_cast<size_t>(m_header.array_size[2]) };
    
    size_t element_size{ GetElementSize() };
    size_t total_bytes{ num_voxels * element_size };
    
    auto data_array{ std::make_unique<float[]>(num_voxels) };
    switch (static_cast<MODE>(m_header.mode))
    {
        case MODE::SIGNED_FLOAT32:
            stream.read(reinterpret_cast<char*>(data_array.get()),
                        static_cast<std::streamsize>(total_bytes));
            if (!stream)
            {
                throw std::runtime_error("Failed to read voxel data from file");
            }
            m_data_array = std::move(data_array);
            break;
        default:
            throw std::runtime_error("Unsupported MODE in LoadDataArray");
    }
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

std::unique_ptr<float[]> MrcFormat::GetDataArray(void)
{
    return std::move(m_data_array);
}

std::array<int, 3> MrcFormat::GetGridSize(void)
{
    std::array<int, 3> grid_size;
    grid_size.at(0) = m_header.grid_size[0];
    grid_size.at(1) = m_header.grid_size[1];
    grid_size.at(2) = m_header.grid_size[2];
    return grid_size;
}

std::array<float, 3> MrcFormat::GetGridSpacing(void)
{
    if (m_header.grid_size[0] == 0 || m_header.grid_size[1] == 0 || m_header.grid_size[2] == 0)
    {
        throw std::runtime_error("GetGridSpacing: grid_size has zero dimension");
    }
    std::array<float, 3> grid_spacing;
    grid_spacing.at(0) = m_header.cell_dimension[0] / static_cast<float>(m_header.grid_size[0]);
    grid_spacing.at(1) = m_header.cell_dimension[1] / static_cast<float>(m_header.grid_size[1]);
    grid_spacing.at(2) = m_header.cell_dimension[2] / static_cast<float>(m_header.grid_size[2]);
    return grid_spacing;
}

std::array<float, 3> MrcFormat::GetOrigin(void)
{
    std::array<float, 3> origin;
    origin.at(0) = m_header.origin[0];
    origin.at(1) = m_header.origin[1];
    origin.at(2) = m_header.origin[2];
    return origin;
}

void MrcFormat::SetGridSize(const std::array<int, 3> & grid_size)
{
    std::memcpy(m_header.array_size, grid_size.data(), sizeof(m_header.array_size));
    std::memcpy(m_header.grid_size, grid_size.data(), sizeof(m_header.grid_size));
}

void MrcFormat::SetGridSpacing(const std::array<float, 3> & grid_spacing)
{
    m_header.cell_dimension[0] = grid_spacing.at(0) * static_cast<float>(m_header.grid_size[0]);
    m_header.cell_dimension[1] = grid_spacing.at(1) * static_cast<float>(m_header.grid_size[1]);
    m_header.cell_dimension[2] = grid_spacing.at(2) * static_cast<float>(m_header.grid_size[2]);
}

void MrcFormat::SetOrigin(const std::array<float, 3> & origin)
{
    std::memcpy(m_header.origin, origin.data(), sizeof(m_header.origin));
}
