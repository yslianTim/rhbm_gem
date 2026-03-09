#define _CRT_SECURE_NO_WARNINGS // To disable deprecation warnings for sscanf
#include "internal/PdbFormat.hpp"
#include <rhbm_gem/data/AtomObject.hpp>
#include <rhbm_gem/utils/ChemicalDataHelper.hpp>
#include <rhbm_gem/data/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/ModelObject.hpp>
#include <rhbm_gem/utils/StringHelper.hpp>
#include <rhbm_gem/data/AtomicModelDataBlock.hpp>
#include <rhbm_gem/utils/Logger.hpp>

#include <fstream>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <stdexcept>

namespace rhbm_gem {

PdbFormat::PdbFormat() :
    m_data_block{ std::make_unique<AtomicModelDataBlock>() }
{

}

PdbFormat::~PdbFormat()
{

}

void PdbFormat::Read(std::istream & stream, const std::string & source_name)
{
    m_data_block = std::make_unique<AtomicModelDataBlock>();
    if (!stream)
    {
        Logger::Log(LogLevel::Error, "Cannot read the file stream: " + source_name);
        throw std::runtime_error("PdbFormat::Read() failed: invalid input stream.");
    }
    LoadAtomSiteData(stream);
}

void PdbFormat::LoadAtomSiteData(std::istream & stream)
{
    std::string line;
    auto model_number{ 1 };
    while (std::getline(stream, line))
    {
        StringHelper::StripCarriageReturn(line);
        char * buffer{ line.data() };
        std::string header_str{ line.substr(0, static_cast<size_t>(header_string_length)) };
        if (header_str == "ENDMDL") model_number++;
        switch (MapToHeaderType(header_str))
        {
        case PDB_HEADER::ATOM:
            ScanAtomEntry(buffer, false, model_number);
            break;
        case PDB_HEADER::HETATM:
            ScanAtomEntry(buffer, true, model_number);
            break;
        default:
            break;
        }
    }
}

void PdbFormat::ScanAtomEntry(char * line, bool is_special, int model_number)
{
    auto atom{ std::make_unique<AtomSite>() };
    std::sscanf(&line[static_cast<int>(ATOM::SERIAL_ID)],    "%d",  &atom->serial_id);
    std::sscanf(&line[static_cast<int>(ATOM::ATOM_NAME)],    "%4c", &atom->atom_name[0]);
    std::sscanf(&line[static_cast<int>(ATOM::INDICATOR)],    "%c",  &atom->indicator);
    std::sscanf(&line[static_cast<int>(ATOM::RESIDUE_NAME)], "%3c", &atom->residue_name[0]);
    std::sscanf(&line[static_cast<int>(ATOM::CHAIN_ID)],     "%c",  &atom->chain_id);
    std::sscanf(&line[static_cast<int>(ATOM::RESIDUE_ID)],   "%d",  &atom->sequence_id);
    std::sscanf(&line[static_cast<int>(ATOM::CODE)],         "%c",  &atom->code);
    std::sscanf(&line[static_cast<int>(ATOM::POSITION_X)],   "%8f", &atom->position_x);
    std::sscanf(&line[static_cast<int>(ATOM::POSITION_Y)],   "%8f", &atom->position_y);
    std::sscanf(&line[static_cast<int>(ATOM::POSITION_Z)],   "%8f", &atom->position_z);
    std::sscanf(&line[static_cast<int>(ATOM::OCCUPANCY)],    "%8f", &atom->occupancy);
    std::sscanf(&line[static_cast<int>(ATOM::TEMPERATURE)],  "%8f", &atom->temperature);
    std::sscanf(&line[static_cast<int>(ATOM::SEGMENT_ID)],   "%4c", &atom->segment_id[0]);
    std::sscanf(&line[static_cast<int>(ATOM::ELEMENT)],      "%2c", &atom->element[0]);
    std::sscanf(&line[static_cast<int>(ATOM::CHARGE)],       "%2c", &atom->charge[0]);

    auto atom_object{ std::make_unique<AtomObject>() };
    auto atom_name{ StringHelper::ConvertCharArrayToString(atom->atom_name) };
    auto element_name{ StringHelper::ConvertCharArrayToString(atom->element) };
    if (element_name == "H") return; // Skip hydrogen atom
    auto indicator{ (atom->indicator == ' ') ? "." : std::string(1, atom->indicator) };
    atom_object->SetElement(element_name);
    atom_object->SetComponentID(atom->segment_id);
    atom_object->SetAtomID(atom_name);
    atom_object->SetResidue(StringHelper::ConvertCharArrayToString(atom->residue_name));
    atom_object->SetIndicator(indicator);
    atom_object->SetSequenceID(atom->sequence_id);
    atom_object->SetSerialID(atom->serial_id);
    atom_object->SetChainID(std::string(1, atom->chain_id));
    atom_object->SetPosition(atom->position_x, atom->position_y, atom->position_z);
    atom_object->SetOccupancy(atom->occupancy);
    atom_object->SetTemperature(atom->temperature);
    atom_object->SetSpecialAtomFlag(is_special);
    m_data_block->AddAtomObject(model_number, std::move(atom_object));
}

AtomicModelDataBlock * PdbFormat::GetDataBlockPtr()
{
    return m_data_block.get();
}

void PdbFormat::Write(const ModelObject & model_object, std::ostream & stream, int par)
{
    if (!stream)
    {
        throw std::runtime_error("PdbFormat::Write() failed: invalid output stream.");
    }

    for (const auto & atom_ptr : model_object.GetAtomList())
    {
        const auto * atom{ atom_ptr.get() };
        if (atom == nullptr)
        {
            continue;
        }

        const auto group_pdb{ atom->GetSpecialAtomFlag() ? "HETATM" : "ATOM" };
        auto atom_name{ atom->GetAtomID() };
        if (atom_name.size() > 4) atom_name = atom_name.substr(0, 4);
        auto component_id{ atom->GetComponentID() };
        auto residue_name{ component_id.substr(0, std::min<size_t>(3, component_id.size())) };
        std::string segment_id{ component_id.substr(0, std::min<size_t>(4, component_id.size())) };
        if (segment_id.empty()) segment_id = "UNK";
        auto element_symbol{ ChemicalDataHelper::GetLabel(atom->GetElement()) };
        if (element_symbol.size() > 2) element_symbol = element_symbol.substr(0, 2);

        const auto indicator{ atom->GetIndicator() };
        const auto alt_loc{
            (indicator.empty() || indicator == ".") ? ' ' : indicator.front()
        };
        const auto chain_id{ atom->GetChainID().empty() ? ' ' : atom->GetChainID().front() };
        const auto sequence_id{ atom->GetSequenceID() < 0 ? 0 : atom->GetSequenceID() };
        const auto position{ atom->GetPosition() };
        const auto occupancy{ atom->GetOccupancy() };
        float b_factor{ atom->GetTemperature() };
        if (atom->GetLocalPotentialEntry() != nullptr)
        {
            b_factor = static_cast<float>(atom->GetLocalPotentialEntry()->GetGausEstimateMDPDE(par));
        }

        char line[128];
        std::snprintf(
            line,
            sizeof(line),
            "%-6s%5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s",
            group_pdb,
            atom->GetSerialID(),
            atom_name.c_str(),
            alt_loc,
            residue_name.c_str(),
            chain_id,
            sequence_id,
            position[0],
            position[1],
            position[2],
            occupancy,
            b_factor,
            segment_id.c_str(),
            element_symbol.c_str());
        stream << line << '\n';
    }
    stream << "END\n";
    if (!stream)
    {
        throw std::runtime_error("PdbFormat::Write() failed while writing output.");
    }
}

PdbFormat::PDB_HEADER PdbFormat::MapToHeaderType(const std::string & name) const
{
    try
    {
        return header_map.at(name);
    }
    catch(const std::exception & except)
    {
        Logger::Log(LogLevel::Error, "Failed to map header type: " + name);
        Logger::Log(LogLevel::Error, "Exception: " + std::string(except.what()));
        return PDB_HEADER::UNK;
    }
}
} // namespace rhbm_gem
