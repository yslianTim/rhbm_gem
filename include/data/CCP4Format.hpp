#pragma once

#include <cstddef>
#include <vector>
#include <istream>
#include <ostream>
#include "MapFileFormatBase.hpp"

class CCP4Format : public MapFileFormatBase
{
    enum HEAD
    {
        SIZE_HEADER = 1024,
        NUM_LABEL   = 10,
        SIZE_LABEL  = 80
    };

    enum AXIS
    {
        COLUMN  = 0,
        ROW     = 1,
        SECTION = 2
    };

    enum class MODE : unsigned int
    {
        SIGNED_INT8     =   0,
        SIGNED_INT16    =   1,
        SIGNED_FLOAT32  =   2,
        COMPLEX_INT16   =   3,
        COMPLEX_FLOAT32 =   4
    };

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
    std::unique_ptr<float[]> m_data_array;

public:
    CCP4Format(void);
    ~CCP4Format() = default;
    void InitHeader(void) override;
    void LoadHeader(std::istream & stream) override;
    void SaveHeader(std::ostream & stream) override;
    void PrintHeader(void) const override;
    void LoadDataArray(std::istream & stream) override;
    void SaveDataArray(const float * data, size_t size, std::ostream & stream) override;
    std::unique_ptr<float[]> GetDataArray(void) override;
    std::array<int, 3> GetGridSize(void) override;
    std::array<float, 3> GetGridSpacing(void) override;
    std::array<float, 3> GetOrigin(void) override;
    void SetGridSize(const std::array<int, 3> & grid_size) override;
    void SetGridSpacing(const std::array<float, 3> & grid_spacing) override;
    void SetOrigin(const std::array<float, 3> & origin) override;

    const CCP4Header & GetHeader(void) const { return m_header; }
    
private:
    size_t GetElementSize(void) const;
    
};
