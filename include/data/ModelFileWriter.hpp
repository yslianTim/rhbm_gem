#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <array>
#include "FileWriterBase.hpp"

class ModelFileFormatBase;
class ModelObject;

class ModelFileWriter : public FileWriterBase
{
    std::string m_file_path;
    std::unique_ptr<ModelFileFormatBase> m_file_format_helper;
    const ModelObject * m_model_object;

public:
    ModelFileWriter(const std::string & filename, const ModelObject * model_object);
    ~ModelFileWriter();
    void Write(void) override;

private:

};