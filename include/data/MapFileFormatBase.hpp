#pragma once

#include <iostream>
#include <cstddef>
#include <memory>
#include <string>
#include <array>

class MapFileFormatBase
{
public:
    virtual ~MapFileFormatBase() = default;
    virtual void InitHeader(void) = 0;
    virtual void LoadHeader(const std::string & filename) = 0;
    virtual void SaveHeader(const std::string & filename) = 0;
    virtual void PrintHeader(void) const = 0;
    virtual void LoadDataArray(const std::string & filename) = 0;
    virtual void SaveDataArray(const std::string & filename) = 0;
    virtual std::unique_ptr<float[]> GetDataArray(void) = 0;
    virtual std::array<int, 3> GetGridSize(void) = 0;
    virtual std::array<float, 3> GetGridSpacing(void) = 0;
    virtual std::array<float, 3> GetOrigin(void) = 0;
    virtual void SetDataArray(size_t array_size, const float * data_array) = 0;
    virtual void SetGridSize(const std::array<int, 3> & grid_size) = 0;
    virtual void SetGridSpacing(const std::array<float, 3> & grid_spacing) = 0;
    virtual void SetOrigin(const std::array<float, 3> & origin) = 0;
};