#include "MrcFormat.hpp"

#include <fstream>
#include <iostream>
#include <cstring>
#include <stdexcept>

MrcFormat::MrcFormat(void)
{
    InitHeader();
}

MrcFormat::~MrcFormat()
{

}

void MrcFormat::InitHeader(void)
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
    m_header.cell_dimension[0] = 1.0f;
    m_header.cell_dimension[1] = 1.0f;
    m_header.cell_dimension[2] = 1.0f;
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
    m_header.extra_size = 0;
    for (int i = 0; i < HEAD::SIZE_WORD; i++) m_header.extra[i] = '0';
    m_header.imodStamp = 0;
    m_header.imodFlags = 0;
    m_header.idtype = 0;
    m_header.lens = 0;
    m_header.nd1 = 0;
    m_header.nd2 = 0;
    m_header.vd1 = 0;
    m_header.vd2 = 0;
    m_header.tiltangles[0] = 0.0f;
    m_header.tiltangles[1] = 0.0f;
    m_header.tiltangles[2] = 0.0f;
    m_header.tiltangles[3] = 0.0f;
    m_header.tiltangles[4] = 0.0f;
    m_header.tiltangles[5] = 0.0f;
    m_header.map[0] = 'M';
    m_header.map[1] = 'A';
    m_header.map[2] = 'P';
    m_header.map[3] = '\0';
    m_header.stamp[0] = '0';
    m_header.stamp[1] = '0';
    m_header.stamp[2] = '0';
    m_header.stamp[3] = '\0';
    m_header.rms = 0.0f;
    m_header.label_size = 0;
    for (int i = 0; i < HEAD::NUM_LABEL; i++)
        for (int j = 0; j <HEAD::SIZE_LABEL; j++)
            m_header.label[i][j] = '0';
}

void MrcFormat::LoadHeader(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadHeader failed!");
    }
    
    infile.read(reinterpret_cast<char*>(&m_header), sizeof(m_header));
    if (!infile)
    {
        throw std::runtime_error("LoadHeader failed!");
    }
}

void MrcFormat::SaveHeader(const std::string & filename)
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

void MrcFormat::PrintHeader(void) const
{
    std::cout << "MRC Header Info:" << std::endl;
    std::cout << "Array Size: " << m_header.array_size[0] << " x "
              << m_header.array_size[1] << " x " << m_header.array_size[2] << std::endl;
    std::cout << "Grid Size: " << m_header.grid_size[0] << " x "
              << m_header.grid_size[1] << " x " << m_header.grid_size[2] << std::endl;
    std::cout << "Cell Dimensions: " << m_header.cell_dimension[0] << ", "
              << m_header.cell_dimension[1] << ", " << m_header.cell_dimension[2] << std::endl;
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

void MrcFormat::LoadDataArray(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        throw std::runtime_error("Cannot open the file: " + filename);
    }
    
    // Skip header (1024 bytes)
    infile.seekg(HEAD::SIZE_HEADER, std::ios::beg);
    
    // Skip if extra header is exist（extra_size > 0)
    if (m_header.extra_size > 0)
    {
        infile.seekg(m_header.extra_size, std::ios::cur);
    }
    
    size_t num_voxels{ static_cast<size_t>(m_header.array_size[0]) *
                       static_cast<size_t>(m_header.array_size[1]) *
                       static_cast<size_t>(m_header.array_size[2]) };
    
    size_t element_size{ GetElementSize() };
    size_t total_bytes{ num_voxels * element_size };
    
    auto blob_buffer{ std::make_unique<char[]>(total_bytes) };
    infile.read(blob_buffer.get(), static_cast<long>(total_bytes));
    if (!infile)
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

void MrcFormat::SaveDataArray(const std::string & filename)
{
    std::fstream file{ filename, std::ios::in | std::ios::out | std::ios::binary };
    if (!file)
    {
        throw std::runtime_error("Cannot open the file: " + filename);
    }

    std::streamoff data_offset{ HEAD::SIZE_HEADER + static_cast<std::streamoff>(m_header.extra_size) };
    file.seekp(data_offset, std::ios::beg);
    if (!file)
    {
        throw std::runtime_error("Failed to seek to data offset");
    }

    const size_t num_voxels{
        static_cast<size_t>(m_header.array_size[0]) *
        static_cast<size_t>(m_header.array_size[1]) *
        static_cast<size_t>(m_header.array_size[2])
    };

    const size_t element_size{ GetElementSize() };
    const size_t total_bytes{ num_voxels * element_size };
    file.write(reinterpret_cast<const char*>(m_data_array.get()),
               static_cast<std::streamsize>(total_bytes));
    if (!file)
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

void MrcFormat::SetDataArray(size_t array_size, const float * data_array)
{
    m_data_array = std::make_unique<float[]>(array_size);
    std::memcpy(m_data_array.get(), data_array, array_size * sizeof(float));
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