#include "PotentialDisplayVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"

PotentialDisplayVisitor::PotentialDisplayVisitor(void)
{

}

PotentialDisplayVisitor::~PotentialDisplayVisitor()
{

}

void PotentialDisplayVisitor::VisitAtomObject(AtomObject * data_object)
{
}

void PotentialDisplayVisitor::VisitModelObject(ModelObject * data_object)
{
    std::cout <<"- PotentialDisplayVisitor::VisitModelObject"<< std::endl;
    
}

void PotentialDisplayVisitor::VisitMapObject(MapObject * data_object)
{
    std::cout <<"This is the map object, do nothing in PotentialDisplayVisitor."<< std::endl;
}

void PotentialDisplayVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- PotentialDisplayVisitor::Analysis" << std::endl;
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