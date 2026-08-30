#pragma once

#include <array>
#include <istream>
#include <memory>
#include <ostream>
#include <string>

namespace rhbm_gem {

class MapObject;

class MrcFormat
{
    enum HEAD
    {
        SIZE_HEADER = 1024,
        NUM_LABEL   = 10,
        SIZE_LABEL  = 80,
        SIZE_WORD   = 56
    };

    static constexpr int kFloat32Mode = 2;

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
        char  map_format_id[4];   // Character string 'MAP' to identify file type
        char  machine_stamp[4];   // Machine stamp encoding byte ordering of data
        float rms;                // RMS deviation of map from mean density
        int   label_size;         // Number of labels being used
        char  label[HEAD::NUM_LABEL][HEAD::SIZE_LABEL];  // 10 80-character text labels
    };
    static_assert(sizeof(MrcHeader) == HEAD::SIZE_HEADER,
                  "MrcHeader size mismatch: check HEAD::SIZE_HEADER");

    MrcHeader m_header;

public:
    MrcFormat();
    ~MrcFormat() = default;
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
