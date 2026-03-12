#include <rhbm_gem/data/io/FileIO.hpp>

#include "internal/io/file/FileFormatBackendFactory.hpp"
#include "internal/io/file/FileFormatRegistry.hpp"
#include "internal/io/file/MapFileFormatBase.hpp"
#include "internal/io/file/ModelFileFormatBase.hpp"
#include <rhbm_gem/data/dispatch/DataObjectDispatch.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomicModelDataBlock.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <fstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace rhbm_gem {

namespace {

const FileFormatRegistry & DefaultFileFormatRegistry()
{
    static const FileFormatRegistry registry{ BuildDefaultFileFormatRegistry() };
    return registry;
}

const FileFormatDescriptor & ResolveDescriptorForRead(const std::filesystem::path & filename)
{
    return DefaultFileFormatRegistry().LookupForRead(FilePathHelper::GetExtension(filename));
}

const FileFormatDescriptor & ResolveDescriptorForWrite(const std::filesystem::path & filename)
{
    return DefaultFileFormatRegistry().LookupForWrite(FilePathHelper::GetExtension(filename));
}

std::unique_ptr<ModelObject> BuildModelObject(AtomicModelDataBlock & data_block)
{
    auto model_number_list{ data_block.GetModelNumberList() };
    if (model_number_list.empty())
    {
        throw std::runtime_error("No atom model found in the input model file.");
    }

    int selected_model_number{ 1 };
    if (data_block.HasModelNumber(selected_model_number) == false)
    {
        selected_model_number = model_number_list.front();
        Logger::Log(LogLevel::Warning,
            "Model 1 not found. Fallback to model "
            + std::to_string(selected_model_number) + ".");
    }

    auto model_object{
        std::make_unique<ModelObject>(data_block.MoveAtomObjectList(selected_model_number))
    };

    std::unordered_set<const AtomObject *> selected_atom_set;
    selected_atom_set.reserve(model_object->GetAtomList().size());
    for (const auto & atom : model_object->GetAtomList())
    {
        selected_atom_set.insert(atom.get());
    }

    auto bond_list{ data_block.MoveBondObjectList() };
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

    model_object->SetPdbID(data_block.GetPdbID());
    model_object->SetEmdID(data_block.GetEmdID());
    model_object->SetResolution(data_block.GetResolution());
    model_object->SetResolutionMethod(data_block.GetResolutionMethod());
    model_object->SetChainIDListMap(data_block.GetChainIDListMap());
    model_object->SetChemicalComponentEntryMap(data_block.GetChemicalComponentEntryMap());
    model_object->SetComponentKeySystem(data_block.MoveComponentKeySystem());
    model_object->SetAtomKeySystem(data_block.MoveAtomKeySystem());
    model_object->SetBondKeySystem(data_block.MoveBondKeySystem());
    model_object->SetBondList(std::move(filtered_bond_list));
    return model_object;
}

std::unique_ptr<ModelObject> ReadModelWithDescriptor(
    const std::filesystem::path & filename,
    const FileFormatDescriptor & descriptor)
{
    if (descriptor.kind != DataObjectKind::Model || !descriptor.model_backend.has_value())
    {
        throw std::runtime_error("Unsupported model file format.");
    }

    auto file_backend{ CreateModelFileFormatBackend(*descriptor.model_backend) };
    std::ifstream file{ filename, std::ios::binary };
    if (!file)
    {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    file_backend->Read(file, filename.string());
    auto * data_block{ file_backend->GetDataBlockPtr() };
    if (data_block == nullptr)
    {
        throw std::runtime_error("Model file backend returned null data block.");
    }
    return BuildModelObject(*data_block);
}

void WriteModelWithDescriptor(
    const std::filesystem::path & filename,
    const FileFormatDescriptor & descriptor,
    const ModelObject & model_object,
    int model_parameter)
{
    if (descriptor.kind != DataObjectKind::Model || !descriptor.model_backend.has_value())
    {
        throw std::runtime_error("Unsupported model file format.");
    }

    auto file_backend{ CreateModelFileFormatBackend(*descriptor.model_backend) };
    std::ofstream outfile{ filename, std::ios::binary };
    if (!outfile)
    {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    file_backend->Write(model_object, outfile, model_parameter);
}

std::unique_ptr<MapObject> ReadMapWithDescriptor(
    const std::filesystem::path & filename,
    const FileFormatDescriptor & descriptor)
{
    if (descriptor.kind != DataObjectKind::Map || !descriptor.map_backend.has_value())
    {
        throw std::runtime_error("Unsupported map file format.");
    }

    auto file_backend{ CreateMapFileFormatBackend(*descriptor.map_backend) };
    std::ifstream file{ filename, std::ios::binary };
    if (!file)
    {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    file_backend->Read(file, filename.string());
    return std::make_unique<MapObject>(
        file_backend->GetGridSize(),
        file_backend->GetGridSpacing(),
        file_backend->GetOrigin(),
        file_backend->GetDataArray());
}

void WriteMapWithDescriptor(
    const std::filesystem::path & filename,
    const FileFormatDescriptor & descriptor,
    const MapObject & map_object)
{
    if (descriptor.kind != DataObjectKind::Map || !descriptor.map_backend.has_value())
    {
        throw std::runtime_error("Unsupported map file format.");
    }

    auto file_backend{ CreateMapFileFormatBackend(*descriptor.map_backend) };
    std::ofstream outfile{ filename, std::ios::binary | std::ios::trunc };
    if (!outfile)
    {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    file_backend->Write(map_object, outfile);
}

} // namespace

std::unique_ptr<DataObjectBase> ReadDataObject(const std::filesystem::path & filename)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForRead(filename) };
        switch (descriptor.kind)
        {
        case DataObjectKind::Model:
            return ReadModelWithDescriptor(filename, descriptor);
        case DataObjectKind::Map:
            return ReadMapWithDescriptor(filename, descriptor);
        }
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to read data object from '" + filename.string() + "': " + ex.what());
    }
    throw std::runtime_error("Failed to read data object from '" + filename.string() + "'.");
}

void WriteDataObject(const std::filesystem::path & filename, const DataObjectBase & data_object)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForWrite(filename) };
        switch (descriptor.kind)
        {
        case DataObjectKind::Model:
            WriteModelWithDescriptor(
                filename,
                descriptor,
                ExpectModelObject(data_object, "WriteDataObject()"),
                0);
            return;
        case DataObjectKind::Map:
            WriteMapWithDescriptor(
                filename,
                descriptor,
                ExpectMapObject(data_object, "WriteDataObject()"));
            return;
        }
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to write data object to '" + filename.string() + "': " + ex.what());
    }
}

std::unique_ptr<ModelObject> ReadModel(const std::filesystem::path & filename)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForRead(filename) };
        return ReadModelWithDescriptor(filename, descriptor);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to read model file '" + filename.string() + "': " + ex.what());
    }
}

void WriteModel(
    const std::filesystem::path & filename,
    const ModelObject & model_object,
    int model_parameter)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForWrite(filename) };
        WriteModelWithDescriptor(filename, descriptor, model_object, model_parameter);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to write model file '" + filename.string() + "': " + ex.what());
    }
}

std::unique_ptr<MapObject> ReadMap(const std::filesystem::path & filename)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForRead(filename) };
        return ReadMapWithDescriptor(filename, descriptor);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to read map file '" + filename.string() + "': " + ex.what());
    }
}

void WriteMap(const std::filesystem::path & filename, const MapObject & map_object)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForWrite(filename) };
        WriteMapWithDescriptor(filename, descriptor, map_object);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to write map file '" + filename.string() + "': " + ex.what());
    }
}

} // namespace rhbm_gem
