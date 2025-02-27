#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <array>

class MapFileFormatBase
{
public:
    virtual ~MapFileFormatBase() = default;
    virtual void LoadHeader(const std::string & filename) = 0;
    virtual void PrintHeader(void) const = 0;
    virtual void LoadDataArray(const std::string & filename) = 0;
    virtual std::unique_ptr<float[]> GetDataArray(void) = 0;
    virtual std::array<int, 3> GetGridSize(void) = 0;
    virtual std::array<float, 3> GetGridSpacing(void) = 0;
    virtual std::array<float, 3> GetOrigin(void) = 0;
};