#include "CCP4Format.hpp"

#include <fstream>
#include <iostream>
#include <cstring>
#include <stdexcept>

CCP4Format::CCP4Format(void)
{
    InitHeader();
}

CCP4Format::~CCP4Format()
{

}

void CCP4Format::InitHeader(void)
{
    std::memset(&m_header, 0, sizeof(m_header));
    m_header.array_size[0] = 1;
    m_header.array_size[1] = 1;
    m_header.array_size[2] = 1;
    m_header.mode = static_cast<int>(MODE::SIGNED_FLOAT32);
    m_header.location_index[0] = 0;
    m_header.location_index[1] = 0;
    m_header.location_index[2] = 0;
    m_header.grid_size[0] = 1;
    m_header.grid_size[1] = 1;
    m_header.grid_size[2] = 1;
    m_header.map_length[0] = 1.0f;
    m_header.map_length[1] = 1.0f;
    m_header.map_length[2] = 1.0f;
    m_header.cell_angle[0] = 90.0f;
    m_header.cell_angle[1] = 90.0f;
    m_header.cell_angle[2] = 90.0f;
    m_header.axis[0] = 1;
    m_header.axis[1] = 2;
    m_header.axis[2] = 3;
    m_header.min_density = 0.0f;
    m_header.max_density = 0.0f;
    m_header.mean_density = 0.0f;
    m_header.space_group = 0;
    m_header.symmetry_table_size = 0;
    m_header.skew_matrix_flag = 0;
    for (int i = 0; i < 9; i++) m_header.skew_matrix[i] = 0.0f;
    for (int i = 0; i < 3; i++) m_header.skew_translation[i] = 0.0f;
    for (int i = 0; i < 15; i++) m_header.extra[i] = 0.0f;
    m_header.map_format_id[0] = 'M';
    m_header.map_format_id[1] = 'A';
    m_header.map_format_id[2] = 'P';
    m_header.map_format_id[3] = '\0';
    m_header.machine_stamp[0] = '0';
    m_header.machine_stamp[1] = '0';
    m_header.machine_stamp[2] = '0';
    m_header.machine_stamp[3] = '\0';
    m_header.rms = 0.0f;
    m_header.label_size = 0;
    for (int i = 0; i < HEAD::NUM_LABEL; i++)
        for (int j = 0; j <HEAD::SIZE_LABEL; j++)
            m_header.label[i][j] = '0';
}

void CCP4Format::LoadHeader(std::istream & stream)
{
    stream.seekg(0, std::ios::beg);
    stream.read(reinterpret_cast<char*>(&m_header), sizeof(m_header));
    if (!stream)
    {
        throw std::runtime_error("LoadHeader failed!");
    }
}

void CCP4Format::SaveHeader(const std::string & filename)
{
    std::ofstream outfile{ filename, std::ios::binary };
    if (!outfile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("SaveHeader failed!");
    }
    
    outfile.write(reinterpret_cast<const char*>(&m_header), sizeof(m_header));
    if (!outfile)
    {
        throw std::runtime_error("SaveHeader failed!");
    }
}

void CCP4Format::PrintHeader(void) const
{
    std::cout << "CCP4 Header Info:" << std::endl;
    std::cout << "Array Size: " << m_header.array_size[0] << " x "
              << m_header.array_size[1] << " x " << m_header.array_size[2] << std::endl;
    std::cout << "Grid Size: " << m_header.grid_size[0] << " x "
              << m_header.grid_size[1] << " x " << m_header.grid_size[2] << std::endl;
    std::cout << "Map Length: " << m_header.map_length[0] << ", "
              << m_header.map_length[1] << ", " << m_header.map_length[2] << std::endl;
}

size_t CCP4Format::GetElementSize(void) const
{
    switch (static_cast<MODE>(m_header.mode))
    {
        case MODE::SIGNED_INT8:
            return 1;
        case MODE::SIGNED_INT16:
            return 2;
        case MODE::SIGNED_FLOAT32:
            return 4;
        case MODE::COMPLEX_INT16:
            return 4;
        case MODE::COMPLEX_FLOAT32:
            return 8;
        default:
            throw std::runtime_error("Unknown mode value");
    }
}

void CCP4Format::LoadDataArray(std::istream & stream)
{
    // Position stream at start of data section
    stream.seekg(HEAD::SIZE_HEADER, std::ios::beg);
    
    size_t num_voxels{ static_cast<size_t>(m_header.array_size[0]) *
                       static_cast<size_t>(m_header.array_size[1]) *
                       static_cast<size_t>(m_header.array_size[2]) };
    
    size_t element_size{ GetElementSize() };
    size_t total_bytes{ num_voxels * element_size };
    
    auto blob_buffer{ std::make_unique<char[]>(total_bytes) };
    stream.read(blob_buffer.get(), static_cast<long>(total_bytes));
    if (!stream)
    {
        throw std::runtime_error("Failed to read blob data from file");
    }
    
    auto data_array{ std::make_unique<float[]>(num_voxels) };
    switch (static_cast<MODE>(m_header.mode))
    {
        case MODE::SIGNED_FLOAT32:
            #ifdef USE_OPENMP
            #pragma omp parallel for num_threads(4)
            #endif
            for (size_t v = 0; v < num_voxels; v++)
            {
                data_array[v] = *reinterpret_cast<const float*>(blob_buffer.get() + v * element_size);
            }
            m_data_array = std::move(data_array);
            break;
        default:
            throw std::runtime_error("Unsupported MODE in LoadDataArray");
    }
}

void CCP4Format::SaveDataArray(const std::string & filename)
{
    if (!m_data_array)
    {
        throw std::runtime_error("SaveDataArray: data array is empty");
    }

    // Open existing file for update -- do NOT truncate header we just saved.
    std::fstream file{ filename, std::ios::in | std::ios::out | std::ios::binary };
    if (!file)
    {
        throw std::runtime_error("Cannot open file for updating: " + filename);
    }

    const std::streamoff data_offset{
        HEAD::SIZE_HEADER + static_cast<std::streamoff>(m_header.symmetry_table_size)
    };

    file.seekp(data_offset, std::ios::beg);
    if (!file)
    {
        throw std::runtime_error("SaveDataArray: failed to seek to data section");
    }

    // Calculate payload size
    const size_t num_voxels{
        static_cast<size_t>(m_header.array_size[0]) *
        static_cast<size_t>(m_header.array_size[1]) *
        static_cast<size_t>(m_header.array_size[2])
    };

    const size_t element_size{ GetElementSize() };
    const size_t total_bytes{ num_voxels * element_size };

    switch (static_cast<MODE>(m_header.mode))
    {
        case MODE::SIGNED_FLOAT32:
        {
            file.write(reinterpret_cast<const char*>(m_data_array.get()),
                       static_cast<std::streamsize>(total_bytes));
            break;
        }
        default:
            throw std::runtime_error("SaveDataArray: unsupported MODE");
    }

    if (!file)
    {
        throw std::runtime_error("SaveDataArray: failed to write voxel data");
    }
}

std::unique_ptr<float[]> CCP4Format::GetDataArray(void)
{
    return std::move(m_data_array);
}

std::array<int, 3> CCP4Format::GetGridSize(void)
{
    std::array<int, 3> grid_size;
    grid_size.at(0) = m_header.grid_size[0];
    grid_size.at(1) = m_header.grid_size[1];
    grid_size.at(2) = m_header.grid_size[2];
    return grid_size;
}

std::array<float, 3> CCP4Format::GetGridSpacing(void)
{
    std::array<float, 3> grid_spacing;
    grid_spacing.at(0) = m_header.map_length[0] / static_cast<float>(m_header.grid_size[0]);
    grid_spacing.at(1) = m_header.map_length[1] / static_cast<float>(m_header.grid_size[1]);
    grid_spacing.at(2) = m_header.map_length[2] / static_cast<float>(m_header.grid_size[2]);
    return grid_spacing;
}

std::array<float, 3> CCP4Format::GetOrigin(void)
{
    std::array<float, 3> origin;
    origin.at(0) = 0.0f;
    origin.at(1) = 0.0f;
    origin.at(2) = 0.0f;
    return origin;
}

void CCP4Format::SetDataArray(size_t array_size, const float * data_array)
{
    m_data_array = std::make_unique<float[]>(array_size);
    std::memcpy(m_data_array.get(), data_array, array_size * sizeof(float));
}

void CCP4Format::SetGridSize(const std::array<int, 3> & grid_size)
{
    std::memcpy(m_header.array_size, grid_size.data(), sizeof(m_header.array_size));
    std::memcpy(m_header.grid_size, grid_size.data(), sizeof(m_header.grid_size));
}

void CCP4Format::SetGridSpacing(const std::array<float, 3> & grid_spacing)
{
    m_header.map_length[0] = grid_spacing.at(0) * static_cast<float>(m_header.grid_size[0]);
    m_header.map_length[1] = grid_spacing.at(1) * static_cast<float>(m_header.grid_size[1]);
    m_header.map_length[2] = grid_spacing.at(2) * static_cast<float>(m_header.grid_size[2]);
}

void CCP4Format::SetOrigin(const std::array<float, 3> & origin)
{
    if (origin.at(0) != 0.0f || origin.at(1) != 0.0f || origin.at(2) != 0.0f)
    {
        throw std::runtime_error("CCP4 format does not support setting origin");
    }
}
