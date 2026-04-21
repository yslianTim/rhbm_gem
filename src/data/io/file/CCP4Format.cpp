#include "CCP4Format.hpp"
#include "MapGridValidation.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include "MapAxisOrderHelper.hpp"

#include <fstream>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace rhbm_gem {

CCP4Format::CCP4Format() {
    InitHeader();
}

std::unique_ptr<MapObject> CCP4Format::ReadMap(
    std::istream& stream, const std::string& source_name) {
    m_data_array.reset();
    InitHeader();
    if (!stream) {
        throw std::runtime_error("CCP4Format::Read() failed: invalid input stream.");
    }

    try {
        LoadHeader(stream);
        PrintHeader();
        LoadDataArray(stream);
    } catch (const std::exception& ex) {
        throw std::runtime_error("CCP4Format::Read() failed for '" + source_name + "': " + ex.what());
    }
    return std::make_unique<MapObject>(GetGridSize(), GetGridSpacing(), GetOrigin(), TakeDataArray());
}

void CCP4Format::WriteMap(const MapObject& map_object, std::ostream& stream) {
    if (!stream) {
        throw std::runtime_error("CCP4Format::Write() failed: invalid output stream.");
    }
    const float* data{map_object.GetMapValueArray()};
    if (data == nullptr) {
        throw std::runtime_error("CCP4Format::Write() failed: map value array is null.");
    }

    InitHeader();
    SetHeader(map_object.GetGridSize(), map_object.GetGridSpacing(), map_object.GetOrigin());
    SaveHeader(stream);
    PrintHeader();
    SaveDataArray(data, map_object.GetMapValueArraySize(), stream);
}

void CCP4Format::InitHeader() {
    std::memset(&m_header, 0, sizeof(m_header));
    std::fill_n(m_header.array_size, 3, 1);
    m_header.mode = static_cast<int>(MODE::SIGNED_FLOAT32);
    std::fill_n(m_header.location_index, 3, 0);
    std::fill_n(m_header.grid_size, 3, 1);
    std::fill_n(m_header.map_length, 3, 1.0f);
    std::fill_n(m_header.cell_angle, 3, 90.0f);
    m_header.axis[0] = 1;
    m_header.axis[1] = 2;
    m_header.axis[2] = 3;
    m_header.min_density = 0.0f;
    m_header.max_density = 0.0f;
    m_header.mean_density = 0.0f;
    m_header.space_group = 0;
    m_header.symmetry_table_size = 0;
    m_header.skew_matrix_flag = 0;
    std::fill_n(m_header.skew_matrix, 9, 0.0f);
    std::fill_n(m_header.skew_translation, 3, 0.0f);
    std::fill_n(m_header.extra, 15, 0.0f);
    m_header.map_format_id[0] = 'M';
    m_header.map_format_id[1] = 'A';
    m_header.map_format_id[2] = 'P';
    m_header.map_format_id[3] = '\0';
    std::fill_n(m_header.machine_stamp, 4, '\0');
    m_header.rms = 0.0f;
    m_header.label_size = 0;
    std::fill_n(&m_header.label[0][0], HEAD::NUM_LABEL * HEAD::SIZE_LABEL, '\0');
}

void CCP4Format::LoadHeader(std::istream& stream) {
    stream.seekg(0, std::ios::beg);
    stream.read(reinterpret_cast<char*>(&m_header), sizeof(m_header));
    if (!stream) {
        throw std::runtime_error("LoadHeader failed!");
    }
}

void CCP4Format::SaveHeader(std::ostream& stream) {
    stream.seekp(0, std::ios::beg);
    stream.write(reinterpret_cast<const char*>(&m_header), sizeof(m_header));
    if (!stream) {
        throw std::runtime_error("SaveHeader failed!");
    }
}

void CCP4Format::PrintHeader() const {
    Logger::Log(LogLevel::Debug,
                "CCP4 Header Information:\n"
                "Array Size: " +
                    std::to_string(m_header.array_size[0]) + " x " + std::to_string(m_header.array_size[1]) + " x " + std::to_string(m_header.array_size[2]) + "\n"
                                                                                                                                                               "Mode: " +
                    std::to_string(m_header.mode) + "\n"
                                                    "Location Index: " +
                    std::to_string(m_header.location_index[0]) + " x " + std::to_string(m_header.location_index[1]) + " x " + std::to_string(m_header.location_index[2]) + "\n"
                                                                                                                                                                           "Grid Size: " +
                    std::to_string(m_header.grid_size[0]) + " x " + std::to_string(m_header.grid_size[1]) + " x " + std::to_string(m_header.grid_size[2]) + "\n"
                                                                                                                                                            "Map Length: " +
                    std::to_string(m_header.map_length[0]) + ", " + std::to_string(m_header.map_length[1]) + ", " + std::to_string(m_header.map_length[2]) + "\n"
                                                                                                                                                             "Cell Angles: " +
                    std::to_string(m_header.cell_angle[0]) + ", " + std::to_string(m_header.cell_angle[1]) + ", " + std::to_string(m_header.cell_angle[2]) + "\n"
                                                                                                                                                             "Axis: " +
                    std::to_string(m_header.axis[0]) + ", " + std::to_string(m_header.axis[1]) + ", " + std::to_string(m_header.axis[2]) + "\n"
                                                                                                                                           "Min Density: " +
                    std::to_string(m_header.min_density) + "\n"
                                                           "Max Density: " +
                    std::to_string(m_header.max_density) + "\n"
                                                           "Mean Density: " +
                    std::to_string(m_header.mean_density) + "\n"
                                                            "Space Group: " +
                    std::to_string(m_header.space_group) + "\n"
                                                           "Symmetry Table Size: " +
                    std::to_string(m_header.symmetry_table_size) + "\n"
                                                                   "Skew Matrix Flag: " +
                    std::to_string(m_header.skew_matrix_flag) + "\n"
                                                                "Skew Matrix: " +
                    std::to_string(m_header.skew_matrix[0]) + ", " + std::to_string(m_header.skew_matrix[1]) + ", " + std::to_string(m_header.skew_matrix[2]) + ", " + std::to_string(m_header.skew_matrix[3]) + ", " + std::to_string(m_header.skew_matrix[4]) + ", " + std::to_string(m_header.skew_matrix[5]) + ", " + std::to_string(m_header.skew_matrix[6]) + ", " + std::to_string(m_header.skew_matrix[7]) + ", " + std::to_string(m_header.skew_matrix[8]) + "\n"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "Skew Translation: " +
                    std::to_string(m_header.skew_translation[0]) + ", " + std::to_string(m_header.skew_translation[1]) + ", " + std::to_string(m_header.skew_translation[2]) + "\n"
                                                                                                                                                                               "Extra: " +
                    std::to_string(m_header.extra[0]) + ", " + std::to_string(m_header.extra[1]) + ", " + std::to_string(m_header.extra[2]) + ", " + std::to_string(m_header.extra[3]) + ", " + std::to_string(m_header.extra[4]) + ", " + std::to_string(m_header.extra[5]) + ", " + std::to_string(m_header.extra[6]) + ", " + std::to_string(m_header.extra[7]) + ", " + std::to_string(m_header.extra[8]) + ", " + std::to_string(m_header.extra[9]) + ", " + std::to_string(m_header.extra[10]) + ", " + std::to_string(m_header.extra[11]) + ", " + std::to_string(m_header.extra[12]) + ", " + std::to_string(m_header.extra[13]) + ", " + std::to_string(m_header.extra[14]) + "\n"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "Map Format ID: " +
                    std::string(m_header.map_format_id) + "\n"
                                                          "Machine Stamp: " +
                    std::string(m_header.machine_stamp) + "\n"
                                                          "RMS: " +
                    std::to_string(m_header.rms) + "\n"
                                                   "Label Size: " +
                    std::to_string(m_header.label_size) + "\n");
}

size_t CCP4Format::GetElementSize() const {
    switch (static_cast<MODE>(m_header.mode)) {
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

void CCP4Format::LoadDataArray(std::istream& stream) {
    // Position stream at start of data section
    stream.seekg(HEAD::SIZE_HEADER, std::ios::beg);
    const std::array<int, 3> array_size{
        m_header.array_size[0],
        m_header.array_size[1],
        m_header.array_size[2]};
    const std::array<int, 3> axis_order{
        m_header.axis[0],
        m_header.axis[1],
        m_header.axis[2]};
    size_t num_voxels{map_io::CountVoxelCount(array_size)};

    size_t element_size{GetElementSize()};
    size_t total_bytes{num_voxels * element_size};

    // Read raw data into temporary buffer first
    auto raw_data{std::make_unique<float[]>(num_voxels)};
    switch (static_cast<MODE>(m_header.mode)) {
    case MODE::SIGNED_FLOAT32:
        stream.read(reinterpret_cast<char*>(raw_data.get()),
                    static_cast<std::streamsize>(total_bytes));
        if (!stream) {
            throw std::runtime_error("Failed to read voxel data from file");
        }
        break;
    default:
        throw std::runtime_error("Unsupported MODE in LoadDataArray");
    }

    auto reordered_array{
        map_io::ReorderToCanonicalXYZ(std::move(raw_data), array_size, axis_order)};

    // Update header to reflect canonical axis order and dimensions
    ReorderedAxisRelatedParameters();

    m_data_array = std::move(reordered_array);
}

void CCP4Format::SaveDataArray(const float* data, size_t size, std::ostream& stream) {
    const std::array<int, 3> array_size{
        m_header.array_size[0],
        m_header.array_size[1],
        m_header.array_size[2]};
    size_t expected_voxels{map_io::CountVoxelCount(array_size)};
    if (size != expected_voxels) {
        throw std::runtime_error("SaveDataArray: voxel count does not match header");
    }

    const std::streamoff data_offset{
        HEAD::SIZE_HEADER + static_cast<std::streamoff>(m_header.symmetry_table_size)};

    stream.seekp(data_offset, std::ios::beg);
    if (!stream) {
        throw std::runtime_error("SaveDataArray: failed to seek to data section");
    }

    const size_t element_size{GetElementSize()};
    const size_t total_bytes{size * element_size};

    switch (static_cast<MODE>(m_header.mode)) {
    case MODE::SIGNED_FLOAT32: {
        stream.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(total_bytes));
        break;
    }
    default:
        throw std::runtime_error("SaveDataArray: unsupported MODE");
    }

    if (!stream) {
        throw std::runtime_error("SaveDataArray: failed to write voxel data");
    }
}

std::unique_ptr<float[]> CCP4Format::TakeDataArray() {
    return std::move(m_data_array);
}

std::array<int, 3> CCP4Format::GetGridSize() const {
    // Return data array size in X, Y, Z order (CCP4Header::array_size)
    std::array<int, 3> grid_size{
        m_header.array_size[0],
        m_header.array_size[1],
        m_header.array_size[2]};
    return grid_size;
}

std::array<float, 3> CCP4Format::GetGridSpacing() const {
    if (m_header.grid_size[0] == 0 || m_header.grid_size[1] == 0 || m_header.grid_size[2] == 0) {
        throw std::runtime_error("GetGridSpacing: grid_size has zero dimension");
    }
    std::array<float, 3> grid_spacing{
        m_header.map_length[0] / static_cast<float>(m_header.grid_size[0]),
        m_header.map_length[1] / static_cast<float>(m_header.grid_size[1]),
        m_header.map_length[2] / static_cast<float>(m_header.grid_size[2])};
    return grid_spacing;
}

std::array<float, 3> CCP4Format::GetOrigin() const {
    auto grid_spacing{GetGridSpacing()};
    std::array<float, 3> origin{
        static_cast<float>(m_header.location_index[0]) * grid_spacing[0],
        static_cast<float>(m_header.location_index[1]) * grid_spacing[1],
        static_cast<float>(m_header.location_index[2]) * grid_spacing[2]};
    return origin;
}

void CCP4Format::SetHeader(const std::array<int, 3>& grid_size,
                           const std::array<float, 3>& grid_spacing,
                           const std::array<float, 3>& origin) {
    NumericValidation::RequireAllPositive(grid_size, "grid_size");
    NumericValidation::RequireAllFinitePositive(grid_spacing, "grid_spacing");

    std::memcpy(m_header.array_size, grid_size.data(), sizeof(m_header.array_size));
    std::memcpy(m_header.grid_size, grid_size.data(), sizeof(m_header.grid_size));

    for (size_t i = 0; i < 3; i++) {
        m_header.map_length[i] = grid_spacing[i] * static_cast<float>(m_header.grid_size[i]);
        m_header.location_index[i] = static_cast<int>(std::lround(origin[i] / grid_spacing[i]));
    }
}

void CCP4Format::ReorderedAxisRelatedParameters() {
    if (m_header.axis[0] == 1 && m_header.axis[1] == 2 && m_header.axis[2] == 3) {
        Logger::Log(LogLevel::Debug,
                    "CCP4Format::ReorderedAxisRelatedParameters : "
                    "Axis already in X->Y->Z order, no need to reorder.");
        return;
    }

    int axis_to_index[3];
    for (int i = 0; i < 3; i++) {
        // axis values are 1-based, convert to 0-based indices
        axis_to_index[m_header.axis[i] - 1] = i;
    }

    int array_size[3]{
        m_header.array_size[0], m_header.array_size[1], m_header.array_size[2]};
    int location_index[3]{
        m_header.location_index[0], m_header.location_index[1], m_header.location_index[2]};
    int grid_size[3]{
        m_header.grid_size[0], m_header.grid_size[1], m_header.grid_size[2]};
    float map_length[3]{
        m_header.map_length[0], m_header.map_length[1], m_header.map_length[2]};
    float cell_angle[3]{
        m_header.cell_angle[0], m_header.cell_angle[1], m_header.cell_angle[2]};

    // Reorder parameters based on current axis mapping
    for (int i = 0; i < 3; i++) {
        m_header.array_size[i] = array_size[axis_to_index[i]];
        m_header.location_index[i] = location_index[axis_to_index[i]];
        m_header.grid_size[i] = grid_size[axis_to_index[i]];
        m_header.map_length[i] = map_length[axis_to_index[i]];
        m_header.cell_angle[i] = cell_angle[axis_to_index[i]];
    }

    // Update axis to canonical X->Y->Z order
    m_header.axis[0] = 1;
    m_header.axis[1] = 2;
    m_header.axis[2] = 3;
}

} // namespace rhbm_gem
