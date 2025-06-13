#pragma once

#include <string>

class AtomObject;
class AtomicModelDataBlock;

class ModelFileFormatBase
{
public:
    virtual ~ModelFileFormatBase() = default;
    virtual void LoadHeader(const std::string & filename) = 0;
    virtual void PrintHeader(void) const = 0;
    virtual void LoadDataArray(const std::string & filename) = 0;
    virtual AtomicModelDataBlock * GetDataBlockPtr(void) = 0;

};