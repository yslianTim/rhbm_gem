#pragma once

#include <memory>
#include <string>

namespace rhbm_gem {

class ModelFileFormatBase;
class AtomObject;
class AtomicModelDataBlock;

class ModelFileReader
{
    std::string m_file_path;
    std::unique_ptr<ModelFileFormatBase> m_file_object;
    bool m_has_loaded_data;

public:
    ModelFileReader(const std::string & filename);
    ~ModelFileReader();
    void Read();

    AtomicModelDataBlock * GetDataBlockPtr();

private:
};

} // namespace rhbm_gem
