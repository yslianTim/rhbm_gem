#pragma once

#include <array>
#include <istream>
#include <memory>
#include <ostream>
#include <string>

namespace rhbm_gem {

class MapObject;

class CCP4Format
{
    enum HEAD
    {
        SIZE_HEADER = 1024,
        NUM_LABEL   = 10,
        SIZE_LABEL  = 80
    };

    static constexpr int kFloat32Mode = 2;

    struct CCP4Header
    {
        int   array_size[3];      // Number of column/row/section in 3D data array
        int   mode;               // Type of data
        int   location_index[3];  // Location of first column/row/section in unit cell
        int   grid_size[3];       // Sampling along X/Y/Z axis of unit cell
        float map_length[3];      // Map lengths along X,Y,Z in angstroms
        float cell_angle[3];      // Unit cell angles in degrees, convention: 90, 90, 90
        int   axis[3];            // Axis corresponding to column/row/section: 1=X, 2=Y, 3=Z
        float min_density;        // Minimum density value
        float max_density;        // Maximum density value
        float mean_density;       // Mean density value
        int   space_group;        // Space group number: 0, 1, 401
        int   symmetry_table_size; // Size of symmetry table in bytes (multiple of 80)
        int   skew_matrix_flag;
        float skew_matrix[9];
        float skew_translation[3];
        float extra[15];          // User-defined metadata
        char  map_format_id[4];   // MRC/CCP4 MAP format identifier
        char  machine_stamp[4];   // machine stamp
        float rms;                // RMS deviation of map from mean density
        int   label_size;         // Number of labels being used
        char  label[HEAD::NUM_LABEL][HEAD::SIZE_LABEL];  // 10 80-character text labels
    };
    static_assert(sizeof(CCP4Header) == HEAD::SIZE_HEADER,
                  "CCP4Header size mismatch: check HEAD::SIZE_HEADER");

    CCP4Header m_header;

public:
    CCP4Format();
    ~CCP4Format() = default;
    std::unique_ptr<MapObject> ReadMap(std::istream & stream, const std::string & source_name);
    void WriteMap(const MapObject & map_object, std::ostream & stream);

private:
    std::array<int, 3> GetGridSize() const;
    std::array<float, 3> GetGridSpacing() const;
    std::array<float, 3> GetOrigin() const;
    void InitHeader();
    void LoadHeader(std::istream & stream);
    void SaveHeader(std::ostream & stream);
    void SetHeader(const std::array<int, 3> & grid_size,
                   const std::array<float, 3> & grid_spacing,
                   const std::array<float, 3> & origin);
};

} // namespace rhbm_gem
