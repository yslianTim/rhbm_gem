#pragma once

#include <memory>
#include <string>
#include <vector>
#include "FileReaderBase.hpp"

class ModelFileFormatBase;
class AtomObject;
class AtomicModelDataBlock;

class ModelFileReader : public FileReaderBase
{
    std::string m_file_path;
    std::unique_ptr<ModelFileFormatBase> m_file_object;

public:
    ModelFileReader(const std::string & filename);
    ~ModelFileReader();
    void Read(void) override;

    AtomicModelDataBlock * GetDataBlockPtr(void);

private:
    void ReadHeader(void);
    void ReadDataArray(void);

};