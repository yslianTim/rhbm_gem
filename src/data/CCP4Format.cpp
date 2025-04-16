#include "CCP4Format.hpp"

#include <fstream>
#include <iostream>
#include <cstring>
#include <stdexcept>

CCP4Format::CCP4Format(void)
{
    std::memset(&m_header, 0, sizeof(m_header));
}

CCP4Format::~CCP4Format()
{

}

void CCP4Format::LoadHeader(const std::string & filename)
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

void CCP4Format::PrintHeader(void) const
{
    std::cout << "CCP4 Header Info:" << std::endl;
    std::cout << "Array Size: " << m_header.array_size[0] << " x "
              << m_header.array_size[1] << " x " << m_header.array_size[2] << std::endl;
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

void CCP4Format::LoadDataArray(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        throw std::runtime_error("Cannot open the file: " + filename);
    }
    
    // Skip header (1024 bytes)
    infile.seekg(HEAD::SIZE_HEADER, std::ios::beg);
    
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
    origin.at(0) = 0.0;
    origin.at(1) = 0.0;
    origin.at(2) = 0.0;
    return origin;
}