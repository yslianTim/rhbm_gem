#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <array>
#include "DataObjectBase.hpp"

class MapObject : public DataObjectBase
{
    std::string m_key_tag;
    int m_voxel_size;
    std::array<int, 3> m_grid_size;
    std::array<float, 3> m_grid_spacing, m_origin, m_map_length;
    std::array<float, 3> m_overflow, m_underflow, m_upper_bound, m_lower_bound;
    std::unique_ptr<float[]> m_map_value_array;
    

public:
    MapObject(const std::array<int, 3> & grid_size,
              const std::array<float, 3> & grid_spacing,
              const std::array<float, 3> & origin,
              std::unique_ptr<float[]> map_value_array);
    ~MapObject();
    MapObject(const MapObject & other);
    std::unique_ptr<DataObjectBase> Clone(void) const override;
    void Display(void) const override;
    void Accept(DataObjectVisitorBase * visitor) override;
    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag(void) const override { return m_key_tag; }

    std::array<float, 3> GetGridSpacing(void) { return m_grid_spacing; }
    std::array<float, 3> GetOrigin(void) { return m_origin; }
    
    std::array<int, 3> GetIndexFromPosition(const std::array<float, 3> & position) const;
    int GetGlobalIndex(int index_x, int index_y, int index_z) const;
    float GetMapValue(int index_x, int index_y, int index_z) const;

private:
    void CheckIndex(int index_x, int index_y, int index_z) const;
    void CheckPosition(const std::array<float, 3> & position) const;

};