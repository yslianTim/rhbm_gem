#pragma once

#include <string>
#include <vector>
#include <memory>
#include "MapFileFormatBase.hpp"

class MrcFormat : public MapFileFormatBase
{
    enum HEAD
    {
        SIZE_HEADER = 1024,
        NUM_LABEL   = 10,
        SIZE_LABEL  = 80,
        SIZE_WORD   = 56
    };

    enum AXIS
    {
        COLUMN  = 0,
        ROW     = 1,
        SECTION = 2
    };

    enum class MODE : unsigned int
    {
        SIGNED_INT8     = 0,
        SIGNED_INT16    = 1,
        SIGNED_FLOAT32  = 2,
        COMPLEX_INT16   = 3,
        COMPLEX_FLOAT32 = 4,
        UNSIGNED_INT16  = 6,
        IEEE754_FLOAT16 = 12,
        TWOPERBYTE_BIT4 = 101
    };

    struct MrcHeader
    {
        int   array_size[3];      // Number of column/row/section in 3D data array
        int   mode;               // Type of data
        int   location_index[3];  // Location of first column/row/section in unit cell
        int   grid_size[3];       // Sampling along X/Y/Z axis of unit cell
        float cell_dimension[3];  // Cell's (full structure) dimensions in angstroms
        float cell_angle[3];      // Cell's (full structure) angles in degrees
        int   axis[3];            // Axis corresponding to column/row/section: 1=X, 2=Y, 3=Z
        float min_density;        // Minimum density value
        float max_density;        // Maximum density value
        float mean_density;       // Mean density value
        int   space_group;        // Space group number: 0, 1, 401
        int   extra_size;         // Size of extended header in bytes
        char  extra[HEAD::SIZE_WORD]; // Extra space used for anything: 0 by default
        int   imodStamp;
        int   imodFlags;
        short idtype;
        short lens;
        short nd1;
        short nd2;
        short vd1;
        short vd2;
        float tiltangles[6];      // 0,1,2 = original:  3,4,5 = current
        float origin[3];          // Phase origin or origin of subvolume
        char  map[4];             // Character string 'MAP' to identify file type
        char  stamp[4];           // Machine stamp encoding byte ordering of data
        float rms;                // RMS deviation of map from mean density
        int   label_size;         // Number of labels being used
        char  label[HEAD::NUM_LABEL][HEAD::SIZE_LABEL];  // 10 80-character text labels
    };

    MrcHeader m_header;
    std::unique_ptr<float[]> m_data_array;

public:
    MrcFormat(void);
    ~MrcFormat();
    void LoadHeader(const std::string & filename) override;
    void PrintHeader(void) const override;
    void LoadDataArray(const std::string & filename) override;
    std::unique_ptr<float[]> GetDataArray(void) override;
    std::array<int, 3> GetGridSize(void) override;
    std::array<float, 3> GetGridSpacing(void) override;
    std::array<float, 3> GetOrigin(void) override;
    const MrcHeader & GetHeader(void) const { return m_header; }
    
private:
    size_t GetElementSize(void) const;
    
};