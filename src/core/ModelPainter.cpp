#include "ModelPainter.hpp"
#include "ModelObject.hpp"
#include "ROOTHelper.hpp"
#include "FilePathHelper.hpp"

#include <TCanvas.h>

ModelPainter::ModelPainter(ModelObject * model) :
    m_model_object{ model }, m_folder_path{ "./" }
{

}

ModelPainter::~ModelPainter()
{

}

void ModelPainter::SetFolder(const std::string & folder_path)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
}

void ModelPainter::Painting(void)
{
    std::cout <<"- ModelPainter::Painting"<<std::endl;
    std::cout <<"  Folder path: "<< m_folder_path << std::endl;
    PaintResidueClassGroupGaus("test1.pdf");
}

void ModelPainter::PaintResidueClassGroupGaus(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueGroupGaus"<<std::endl;
    std::cout <<"  Output file: "<< file_path << std::endl;
    auto canvas{ ROOTHelper::CreateCanvas("test","") };

    ROOTHelper::SetCanvasStyle(canvas.get());
    ROOTHelper::SetPadInCanvas(canvas.get(), 1000, 600, 2, 2);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
}