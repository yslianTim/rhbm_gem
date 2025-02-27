#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <array>
#include <vector>
#include "FileReaderBase.hpp"

class ModelFileFormatBase;
class AtomObject;

class ModelFileReader : public FileReaderBase
{
    std::string m_file_path;
    std::unique_ptr<ModelFileFormatBase> m_file_format_helper;

public:
    ModelFileReader(const std::string & filename);
    ~ModelFileReader();
    void Read(void) override;

    std::vector<std::unique_ptr<AtomObject>> GetAtomObjectList(void);
    std::string GetPdbID(void) const;
    std::string GetEmdID(void) const;

private:
    void ReadHeader(void);
    void ReadDataArray(void);

};