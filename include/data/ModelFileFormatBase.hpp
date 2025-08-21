#pragma once

#include <string>
#include <ostream>

class AtomObject;
class AtomicModelDataBlock;
class ModelObject;

class ModelFileFormatBase
{
public:
    virtual ~ModelFileFormatBase() = default;
    virtual void LoadHeader(const std::string & filename) = 0;
    virtual void PrintHeader(void) const = 0;
    virtual void LoadDataArray(const std::string & filename) = 0;
    virtual void SaveHeader(const ModelObject * model_object, std::ostream & stream) = 0;
    virtual void SaveDataArray(const ModelObject * model_object, std::ostream & stream, int par) = 0;
    virtual AtomicModelDataBlock * GetDataBlockPtr(void) = 0;

};
