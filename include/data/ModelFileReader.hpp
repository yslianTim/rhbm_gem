#pragma once

#include <memory>
#include <string>
#include <vector>
#include "FileReaderBase.hpp"

namespace rhbm_gem {

class ModelFileFormatBase;
class AtomObject;
class AtomicModelDataBlock;

class ModelFileReader : public FileReaderBase
{
    bool m_successfully_read_file;
    std::string m_file_path;
    std::unique_ptr<ModelFileFormatBase> m_file_object;

public:
    ModelFileReader(const std::string & filename);
    ~ModelFileReader();
    void Read() override;
    bool IsSuccessfullyRead() const override { return m_successfully_read_file; }

    AtomicModelDataBlock * GetDataBlockPtr();

private:
    void ReadHeader();
    void ReadDataArray();

};

} // namespace rhbm_gem
