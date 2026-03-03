#pragma once

#include <memory>
#include <string>
#include "FileWriterBase.hpp"

namespace rhbm_gem {

class ModelFileFormatBase;
class ModelObject;

class ModelFileWriter : public FileWriterBase
{
    std::string m_file_path;
    std::unique_ptr<ModelFileFormatBase> m_file_object;
    const ModelObject * m_model_object;
    int m_model_par;

public:
    ModelFileWriter(const std::string & filename, const ModelObject * model_object, int model_par = 0);
    ~ModelFileWriter();
    void Write(void) override;

};

} // namespace rhbm_gem
