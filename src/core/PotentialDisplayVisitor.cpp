#include "PotentialDisplayVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "DataObjectBase.hpp"
#include "PainterBase.hpp"
#include "AtomPainter.hpp" 
#include "AtomSelector.hpp"
#include "ModelPainter.hpp"
#include "ComparisonPainter.hpp"
#include "DemoPainter.hpp"
#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"

#include <iostream>
#include <memory>

PotentialDisplayVisitor::PotentialDisplayVisitor(AtomSelector * atom_selector) :
    m_atom_selector{ atom_selector }
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
    BuildReferenceModelObjectList(data_manager, "amber95", m_additional_model_object_list_map["amber95"]);
/* 
    std::vector<std::string> element_list{"Ap", "Cp", "Nn", "On"};
    std::vector<std::string> charge_list{"1", "3", "5", "7"};
    for (auto & element : element_list)
    {
        for (auto & charge : charge_list)
        {
            auto model_class{ element + charge };
            BuildReferenceModelObjectList(data_manager, model_class, m_additional_model_object_list_map[model_class]);
        }
    }
*/
    RunAtomPainter(dynamic_cast<ModelObject *>(m_model_object_list.at(0)));
    RunModelPainter(data_manager);
    RunComparisonPainter(data_manager);
    RunDemoPainter(data_manager);
}

void PotentialDisplayVisitor::RunAtomPainter(ModelObject * model_object)
{
    std::cout <<"- PotentialDisplayVisitor::RunAtomPainter" << std::endl;
    std::unique_ptr<PainterBase> painter{ std::make_unique<AtomPainter>() };
    painter->SetFolder(m_folder_path);
    for (auto & atom : model_object->GetComponentsList())
    {
        if (atom->GetAtomicChargeEntry() == nullptr) continue;
        if (atom->GetRemoteness() != Remoteness::NONE && 
            atom->GetRemoteness() != Remoteness::ALPHA) continue;
        if (atom->GetSpecialAtomFlag() == true) continue;
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

    for (auto & [model_class, model_object_list] : m_additional_model_object_list_map)
    {
        for (auto & model : model_object_list)
        {
            painter->AddReferenceDataObject(model, model_class);
        }
    }
    painter->Painting();
}

void PotentialDisplayVisitor::RunDemoPainter(DataObjectManager * data_manager)
{
    std::cout <<"- PotentialDisplayVisitor::RunDemoPainter" << std::endl;
    if (data_manager == nullptr) return;
    std::unique_ptr<PainterBase> painter{ std::make_unique<DemoPainter>() };
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
