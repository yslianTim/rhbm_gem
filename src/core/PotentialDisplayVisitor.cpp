#include "PotentialDisplayVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
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

PotentialDisplayVisitor::PotentialDisplayVisitor(
    AtomSelector * atom_selector, const PotentialDisplayCommand::Options & options) :
    m_painter_choice{ options.painter_choice },
    m_folder_path{ options.folder_path },
    m_model_key_tag_list{ options.model_key_tag_list },
    m_atom_selector{ atom_selector }
{

}

PotentialDisplayVisitor::~PotentialDisplayVisitor()
{
    BuildOrderedModelObjectList();
    BuildOrderedRefModelObjectListMap();
    std::unique_ptr<PainterBase> painter{ nullptr };
    switch (m_painter_choice)
    {
        case 0:
            painter = std::make_unique<AtomPainter>();
            for (auto * model_object : m_ordered_model_object_list)
            {
                for (auto & atom : model_object->GetComponentsList())
                {
                    if (atom->GetSelectedFlag() == false) continue;
                    painter->AddDataObject(atom.get());
                }
            }
            m_atom_selector->Print();
            break;
        case 1:
            painter = std::make_unique<ModelPainter>();
            for (auto * model_object : m_ordered_model_object_list)
            {
                painter->AddDataObject(model_object);
            }
            for (auto & [class_key, obj_list] : m_ordered_ref_model_object_list_map)
            {
                for (auto * obj : obj_list)
                {
                    painter->AddReferenceDataObject(obj, class_key);
                }
            }
            break;
        case 2:
            painter = std::make_unique<ComparisonPainter>();
            for (auto * model_object : m_ordered_model_object_list)
            {
                painter->AddDataObject(model_object);
            }
            for (auto & [class_key, obj_list] : m_ordered_ref_model_object_list_map)
            {
                for (auto * obj : obj_list)
                {
                    painter->AddReferenceDataObject(obj, class_key);
                }
            }
            break;
        case 3:
            painter = std::make_unique<DemoPainter>();
            for (auto * model_object : m_ordered_model_object_list)
            {
                painter->AddDataObject(model_object);
            }
            for (auto & [class_key, obj_list] : m_ordered_ref_model_object_list_map)
            {
                for (auto * obj : obj_list)
                {
                    painter->AddReferenceDataObject(obj, class_key);
                }
            }
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
    if (painter)
    {
        painter->SetFolder(m_folder_path);
        painter->Painting();
    }
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
    auto key_tag{ data_object->GetKeyTag() };
    for (auto & tag : m_model_key_tag_list)
    {
        if (key_tag == tag)
        {
            m_model_object_list.emplace_back(data_object);
            return;
        }
    }
    for (auto & [class_key, tag_list] : m_ref_model_key_tag_list_map)
    {
        for (auto & tag : tag_list)
        {
            if (key_tag == tag)
            {
                m_ref_model_object_list_map[class_key].emplace_back(data_object);
                return;
            }
        }
    }
}

void PotentialDisplayVisitor::BuildOrderedModelObjectList(void)
{
    m_ordered_model_object_list.clear();
    for (const auto & key : m_model_key_tag_list)
    {
        for (auto * model_object : m_model_object_list)
        {
            if (model_object->GetKeyTag() == key)
            {
                m_ordered_model_object_list.emplace_back(model_object);
                break;
            }
        }
    }
}

void PotentialDisplayVisitor::BuildOrderedRefModelObjectListMap(void)
{
    m_ordered_ref_model_object_list_map.clear();
    for (auto & [class_key, ordered_key_list] : m_ref_model_key_tag_list_map)
    {
        auto it{ m_ref_model_object_list_map.find(class_key) };
        if (it == m_ref_model_object_list_map.end()) continue;
        auto & model_object_list{ it->second };
        std::vector<ModelObject *> ordered_model_object_list;
        for (const auto & key : ordered_key_list)
        {
            for (auto * model_object : model_object_list)
            {
                if (model_object->GetKeyTag() == key)
                {
                    ordered_model_object_list.emplace_back(model_object);
                    break;
                }
            }
        }

        m_ordered_ref_model_object_list_map[class_key] = std::move(ordered_model_object_list);
    }
}
