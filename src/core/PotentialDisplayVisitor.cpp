#include "PotentialDisplayVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "DataObjectBase.hpp"
#include "PainterBase.hpp"
#include "AtomPainter.hpp" 
#include "ModelPainter.hpp"
#include "ComparisonPainter.hpp"
#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"

#include <memory>

PotentialDisplayVisitor::PotentialDisplayVisitor(void)
{

}

PotentialDisplayVisitor::~PotentialDisplayVisitor()
{

}

void PotentialDisplayVisitor::VisitAtomObject(AtomObject * data_object)
{
    if (data_object == nullptr) return;
}

void PotentialDisplayVisitor::VisitModelObject(ModelObject * data_object)
{
    if (data_object == nullptr) return;
}

void PotentialDisplayVisitor::VisitMapObject(MapObject * data_object)
{
    if (data_object == nullptr) return;
}

void PotentialDisplayVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- PotentialDisplayVisitor::Analysis" << std::endl;
    BuildModelObjectList(data_manager, m_model_object_list);
    BuildReferenceModelObjectList(data_manager, "no_charge", m_sim_no_charge_model_object_list);
    BuildReferenceModelObjectList(data_manager, "with_charge", m_sim_with_charge_model_object_list);
    //RunAtomPainter(dynamic_cast<ModelObject *>(m_model_object_list.at(0)));
    RunModelPainter(data_manager);
    RunComparisonPainter(data_manager);
}

void PotentialDisplayVisitor::RunAtomPainter(ModelObject * model_object)
{
    std::cout <<"- PotentialDisplayVisitor::RunAtomPainter" << std::endl;
    std::unique_ptr<PainterBase> painter{ std::make_unique<AtomPainter>() };
    painter->SetFolder(m_folder_path);
    for (auto & atom : model_object->GetComponentsList())
    {
        if (atom->GetSelectedFlag() == false) continue;
        if (atom->GetElement() != Element::CARBON) continue;
        if (atom->GetRemoteness() != Remoteness::ALPHA) continue;
        if (atom->GetResidue() != Residue::ALA) continue;
        painter->AddDataObject(atom.get());
    }
    painter->Painting();
}

void PotentialDisplayVisitor::RunModelPainter(DataObjectManager * data_manager)
{
    std::cout <<"- PotentialDisplayVisitor::RunModelPainter" << std::endl;
    if (data_manager == nullptr) return;
    std::unique_ptr<PainterBase> painter{ std::make_unique<ModelPainter>() };
    painter->SetFolder(m_folder_path);
    for (auto model_object : m_model_object_list)
    {
        painter->AddDataObject(model_object);
    }
    
    for (auto model_object : m_sim_no_charge_model_object_list)
    {
        painter->AddReferenceDataObject(model_object, "no_charge");
    }

    for (auto model_object : m_sim_with_charge_model_object_list)
    {
        painter->AddReferenceDataObject(model_object, "with_charge");
    }
    painter->Painting();
}

void PotentialDisplayVisitor::RunComparisonPainter(DataObjectManager * data_manager)
{
    std::cout <<"- PotentialDisplayVisitor::RunComparisonPainter" << std::endl;
    if (data_manager == nullptr) return;
    std::unique_ptr<PainterBase> painter{ std::make_unique<ComparisonPainter>() };
    painter->SetFolder(m_folder_path);

    for (auto model_object : m_model_object_list)
    {
        painter->AddDataObject(model_object);
    }

    for (auto model_object : m_sim_no_charge_model_object_list)
    {
        painter->AddReferenceDataObject(model_object, "no_charge");
    }

    for (auto model_object : m_sim_with_charge_model_object_list)
    {
        painter->AddReferenceDataObject(model_object, "with_charge");
    }
    painter->Painting();
}

void PotentialDisplayVisitor::BuildModelObjectList(
    DataObjectManager * data_manager, std::vector<DataObjectBase *> & model_object_list)
{
    model_object_list.clear();
    for (auto & key_tag : m_model_key_tag_list)
    {
        try
        {
            model_object_list.emplace_back(data_manager->GetDataObjectRef(key_tag).get());
        }
        catch(const std::exception & e)
        {
            std::cerr << e.what() << std::endl;
        }
    }
}

void PotentialDisplayVisitor::BuildReferenceModelObjectList(
    DataObjectManager * data_manager,
    const std::string & type_key, std::vector<DataObjectBase *> & model_object_list)
{
    model_object_list.clear();
    if (m_ref_model_key_tag_list_map.find(type_key) == m_ref_model_key_tag_list_map.end())
    {
        std::cout <<" Cannot find simulation type : ["<< type_key <<"] in the list."<< std::endl;
        return;
    }
    for (auto & key_tag : m_ref_model_key_tag_list_map.at(type_key))
    {
        try
        {
            model_object_list.emplace_back(data_manager->GetDataObjectRef(key_tag).get());
        }
        catch(const std::exception & e)
        {
            std::cerr << e.what() << std::endl;
        }
    }
}
