#pragma once

#include "PainterBase.hpp"

class ModelObject;

class ModelPainter : public PainterBase
{
    ModelObject * m_model_object;
    std::string m_folder_path;

public:
    ModelPainter(void) = delete;
    explicit ModelPainter(ModelObject * model);
    ~ModelPainter();
    void SetFolder(std::string folder_path) override;
    void Painting(void) override;
};