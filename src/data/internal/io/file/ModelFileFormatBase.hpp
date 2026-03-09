#pragma once

#include <istream>
#include <ostream>
#include <string>

namespace rhbm_gem {

class AtomObject;
class AtomicModelDataBlock;
class ModelObject;

class ModelFileFormatBase
{
public:
    virtual ~ModelFileFormatBase() = default;
    virtual void Read(std::istream & stream, const std::string & source_name) = 0;
    virtual void Write(const ModelObject & model_object, std::ostream & stream, int par) = 0;
    virtual AtomicModelDataBlock * GetDataBlockPtr() = 0;
};

} // namespace rhbm_gem
