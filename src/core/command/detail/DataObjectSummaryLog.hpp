#pragma once

#include <iomanip>
#include <sstream>
#include <string>

#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

namespace rhbm_gem::command_detail {

inline void LogModelSummary(const ModelObject & model_object)
{
    Logger::Log(LogLevel::Info, "ModelObject Summary: " + model_object.GetKeyTag());
    Logger::Log(LogLevel::Info,
        " - PDB ID = " + model_object.GetPdbID() + "\n"
        " - EMD ID = " + model_object.GetEmdID() + "\n"
        " - Map Resolution = "
            + StringHelper::ToStringWithPrecision(model_object.GetResolution(), 2)
            + " A (" + model_object.GetResolutionMethod() + ")\n"
        " - #Unique Entities = " + std::to_string(model_object.GetChainIDListMap().size()) + "\n"
        " - #Atoms = " + std::to_string(model_object.GetNumberOfAtom()) + "\n"
        " - #Bonds = " + std::to_string(model_object.GetNumberOfBond()) + "\n"
        " - #Unique Components = "
            + std::to_string(model_object.GetChemicalComponentEntryMap().size()) + "\n");

    for (const auto & [component_key, entry] : model_object.GetChemicalComponentEntryMap())
    {
        (void)component_key;
        Logger::Log(LogLevel::Debug,
            "   - Component ID = " + entry->GetComponentId() + " "
            " (#Atoms = " + std::to_string(entry->AtomEntries().size()) +
            ", #Bonds = " + std::to_string(entry->BondEntries().size()) + ")");
    }
}

inline void LogMapSummary(const MapObject & map_object)
{
    std::ostringstream oss;
    oss << "MapObject Summary:\n";
    oss << " o=====================================================o\n";
    oss << " |  Map Object  |   X-axis   |   Y-axis   |   Z-axis   |\n";
    oss << " o=====================================================o\n";
    oss << " | Grid size    | ";
    oss << std::setw(10) << map_object.GetGridSize().at(0) << " | "
        << std::setw(10) << map_object.GetGridSize().at(1) << " | "
        << std::setw(10) << map_object.GetGridSize().at(2) << " |\n";
    oss << " | Grid Spacing | ";
    oss << std::setw(10) << map_object.GetGridSpacing().at(0) << " | "
        << std::setw(10) << map_object.GetGridSpacing().at(1) << " | "
        << std::setw(10) << map_object.GetGridSpacing().at(2) << " |\n";
    oss << " | Origin (A)   | ";
    oss << std::setw(10) << map_object.GetOrigin().at(0) << " | "
        << std::setw(10) << map_object.GetOrigin().at(1) << " | "
        << std::setw(10) << map_object.GetOrigin().at(2) << " |\n";
    oss << " | Map Length(A)| ";
    oss << std::setw(10)
        << static_cast<float>(map_object.GetGridSize().at(0)) * map_object.GetGridSpacing().at(0)
        << " | "
        << std::setw(10)
        << static_cast<float>(map_object.GetGridSize().at(1)) * map_object.GetGridSpacing().at(1)
        << " | "
        << std::setw(10)
        << static_cast<float>(map_object.GetGridSize().at(2)) * map_object.GetGridSpacing().at(2)
        << " |\n";
    oss << " |-----------------------------------------------------|\n";
    oss << " | Map value min  | " << std::setw(34) << map_object.GetMapValueMin() << " |\n";
    oss << " | Map value max  | " << std::setw(34) << map_object.GetMapValueMax() << " |\n";
    oss << " | Map value mean | " << std::setw(34) << map_object.GetMapValueMean() << " |\n";
    oss << " | Map value s.d. | " << std::setw(34) << map_object.GetMapValueSD() << " |\n";
    oss << " o=====================================================o\n";
    Logger::Log(LogLevel::Info, oss.str());
}

} // namespace rhbm_gem::command_detail
