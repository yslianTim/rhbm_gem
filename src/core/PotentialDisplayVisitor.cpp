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
    bool selected_flag
    {
        m_atom_selector->GetSelectionFlag(
            data_object->GetChainID(),
            data_object->GetResidue(),
            data_object->GetElement(),
            data_object->GetRemoteness()
        )
    };
    data_object->SetSelectedFlag(selected_flag);
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
    if (data_manager == nullptr) return;

    std::unique_ptr<PainterBase> painter{ nullptr };
    
    switch (m_painter_choice)
    {
        case 0:
            painter = std::make_unique<AtomPainter>();
            AddAtomObjectToPainter(data_manager, painter.get());
            m_atom_selector->Print();
            break;
        case 1:
            painter = std::make_unique<ModelPainter>();
            AddModelObjectToPainter(data_manager, painter.get());
            AddReferenceModelObjectToPainter(data_manager, painter.get());
            break;
        case 2:
            painter = std::make_unique<ComparisonPainter>();
            AddModelObjectToPainter(data_manager, painter.get());
            AddReferenceModelObjectToPainter(data_manager, painter.get());
            break;
        case 3:
            painter = std::make_unique<DemoPainter>();
            AddModelObjectToPainter(data_manager, painter.get());
            AddReferenceModelObjectToPainter(data_manager, painter.get());
            break;
        default:
            std::cout <<"[Warning] Invalid painter choice input : ["<< m_painter_choice <<"]"<< std::endl;
            std::cout <<" * Available Painter Choices : "<< std::endl;
            std::cout <<"   - [0] AtomPainter : "<< std::endl;
            std::cout <<"   - [1] ModelPainter : "<< std::endl;
            std::cout <<"   - [2] ComparisonPainter : "<< std::endl;
            std::cout <<"   - [3] DemoPainter : "<< std::endl;
            break;
    }

    painter->SetFolder(m_folder_path);
    painter->Painting();
}

void PotentialDisplayVisitor::AddAtomObjectToPainter(
    DataObjectManager * data_manager, PainterBase * painter)
{
    for (auto & key_tag : m_model_key_tag_list)
    {
        auto & data_object{ data_manager->GetDataObjectRef(key_tag) };
        auto model_object{ dynamic_cast<ModelObject *>(data_object.get()) };
        model_object->Accept(this);
        for (auto & atom : model_object->GetComponentsList())
        {
            if (atom->GetSelectedFlag() == false) continue;
            painter->AddDataObject(atom.get());
        }
    }
}

void PotentialDisplayVisitor::AddModelObjectToPainter(
    DataObjectManager * data_manager, PainterBase * painter)
{
    for (auto & key_tag : m_model_key_tag_list)
    {
        auto & data_object{ data_manager->GetDataObjectRef(key_tag) };
        painter->AddDataObject(data_object.get());
    }
}

void PotentialDisplayVisitor::AddReferenceModelObjectToPainter(
    DataObjectManager * data_manager, PainterBase * painter)
{
    for (auto & [class_key, key_tag_list] : m_ref_model_key_tag_list_map)
    {
        for (auto & key_tag : key_tag_list)
        {
            auto & data_object{ data_manager->GetDataObjectRef(key_tag) };
            painter->AddReferenceDataObject(data_object.get(), class_key);
        }
    }
}
