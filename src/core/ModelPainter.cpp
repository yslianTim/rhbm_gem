#include "ModelPainter.hpp"
#include "ModelObject.hpp"

ModelPainter::ModelPainter(ModelObject * model) :
    m_model_object{ model }, m_folder_path{ "./" }
{

}

ModelPainter::~ModelPainter()
{

}

void ModelPainter::SetFolder(std::string folder_path)
{
    m_folder_path = folder_path;
}

void ModelPainter::Painting(void)
{
    std::cout <<"- ModelPainter::Painting"<<std::endl;
}