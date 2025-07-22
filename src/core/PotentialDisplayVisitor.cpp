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
#include "Logger.hpp"

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

void PotentialDisplayVisitor::VisitDataObjectManager(DataObjectManager * data_manager)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayVisitor::VisitDataObjectManager() called");
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
            Logger::Log(LogLevel::Warning,
                        "Invalid painter choice input: [" + std::to_string(m_painter_choice) + "]");
            Logger::Log(LogLevel::Warning,
                        "Available Painter Choices:\n"
                        "  [0] AtomPainter\n"
                        "  [1] ModelPainter\n"
                        "  [2] ComparisonPainter\n"
                        "  [3] DemoPainter");
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
        auto data_object{ data_manager->GetDataObjectRef(key_tag) };
        auto model_object{ dynamic_cast<ModelObject *>(data_object) };
        if (model_object == nullptr)
        {
            throw std::runtime_error(
                "PotentialDisplayVisitor::Visit(): invalid model object");
        }
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
        auto data_object{ data_manager->GetDataObjectRef(key_tag) };
        painter->AddDataObject(data_object);
    }
}

void PotentialDisplayVisitor::AddReferenceModelObjectToPainter(
    DataObjectManager * data_manager, PainterBase * painter)
{
    for (auto & [class_key, key_tag_list] : m_ref_model_key_tag_list_map)
    {
        for (auto & key_tag : key_tag_list)
        {
            auto data_object{ data_manager->GetDataObjectRef(key_tag) };
            painter->AddReferenceDataObject(data_object, class_key);
        }
    }
}
