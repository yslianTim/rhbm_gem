#include "PotentialComparisonVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "PainterBase.hpp"
#include "ModelPainter.hpp"

#include <memory>

PotentialComparisonVisitor::PotentialComparisonVisitor(void)
{

}

PotentialComparisonVisitor::~PotentialComparisonVisitor()
{

}

void PotentialComparisonVisitor::VisitAtomObject(AtomObject * data_object)
{
}

void PotentialComparisonVisitor::VisitModelObject(ModelObject * data_object)
{
    std::cout <<"- PotentialComparisonVisitor::VisitModelObject"<< std::endl;
    //auto model_painter{ std::make_unique<ModelPainter>(data_object) };
    //model_painter->SetFolder(m_folder_path);
    //model_painter->Painting();
}

void PotentialComparisonVisitor::VisitMapObject(MapObject * data_object)
{
    std::cout <<"This is the map object, do nothing in PotentialComparisonVisitor."<< std::endl;
}

void PotentialComparisonVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- PotentialComparisonVisitor::Analysis" << std::endl;
    try
    {
        const auto & model_object{ data_manager->GetDataObjectRef(m_model_key_tag) };
        model_object->Accept(this);
    }
    catch(const std::exception & e)
    {
        std::cerr << e.what() << std::endl;
    }
}