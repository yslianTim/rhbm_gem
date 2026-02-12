#include "FileProcessFactoryBase.hpp"
#include "ModelFileReader.hpp"
#include "ModelFileWriter.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "AtomicModelDataBlock.hpp"
#include "Logger.hpp"

#include <stdexcept>
#include <unordered_set>

std::unique_ptr<DataObjectBase> ModelObjectFactory::CreateDataObject(const std::string & filename)
{
    Logger::Log(LogLevel::Debug, "ModelObjectFactory::CreateDataObject() called");
    auto file_reader{ std::make_unique<ModelFileReader>(filename) };
    file_reader->Read();
    if (file_reader->IsSuccessfullyRead() == false) return nullptr;
    auto data_block{ file_reader->GetDataBlockPtr() };
    auto model_number_list{ data_block->GetModelNumberList() };
    if (model_number_list.empty())
    {
        throw std::runtime_error("No atom model found in the input model file.");
    }

    int selected_model_number{ 1 };
    if (data_block->HasModelNumber(selected_model_number) == false)
    {
        selected_model_number = model_number_list.front();
        Logger::Log(LogLevel::Warning,
            "Model 1 not found. Fallback to model "
            + std::to_string(selected_model_number) + ".");
    }

    auto model_object{
        std::make_unique<ModelObject>(data_block->MoveAtomObjectList(selected_model_number))
    };

    std::unordered_set<const AtomObject *> selected_atom_set;
    selected_atom_set.reserve(model_object->GetAtomList().size());
    for (const auto & atom : model_object->GetAtomList())
    {
        selected_atom_set.insert(atom.get());
    }

    auto bond_list{ data_block->MoveBondObjectList() };
    std::vector<std::unique_ptr<BondObject>> filtered_bond_list;
    filtered_bond_list.reserve(bond_list.size());
    for (auto & bond : bond_list)
    {
        if (bond == nullptr) continue;
        auto atom_1{ bond->GetAtomObject1() };
        auto atom_2{ bond->GetAtomObject2() };
        if (selected_atom_set.find(atom_1) == selected_atom_set.end()) continue;
        if (selected_atom_set.find(atom_2) == selected_atom_set.end()) continue;
        filtered_bond_list.emplace_back(std::move(bond));
    }

    model_object->SetPdbID(data_block->GetPdbID());
    model_object->SetEmdID(data_block->GetEmdID());
    model_object->SetResolution(data_block->GetResolution());
    model_object->SetResolutionMethod(data_block->GetResolutionMethod());
    model_object->SetChainIDListMap(data_block->GetChainIDListMap());
    model_object->SetChemicalComponentEntryMap(data_block->GetChemicalComponentEntryMap());
    model_object->SetComponentKeySystem(data_block->MoveComponentKeySystem());
    model_object->SetAtomKeySystem(data_block->MoveAtomKeySystem());
    model_object->SetBondList(std::move(filtered_bond_list));
    return model_object;
}

void ModelObjectFactory::OutputDataObject(const std::string & filename, DataObjectBase * data_object)
{
    Logger::Log(LogLevel::Debug, "ModelObjectFactory::OutputDataObject() called");
    auto model_object{ dynamic_cast<ModelObject *>(data_object) };
    if (model_object == nullptr)
    {
        throw std::runtime_error(
            "ModelObjectFactory::OutputDataObject(): invalid data_object type");
    }
    auto file_writer{ std::make_unique<ModelFileWriter>(filename, model_object) };
    file_writer->Write();
}
