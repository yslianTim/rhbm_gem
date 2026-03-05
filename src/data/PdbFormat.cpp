#define _CRT_SECURE_NO_WARNINGS // To disable deprecation warnings for sscanf
#include "PdbFormat.hpp"
#include "AtomObject.hpp"
#include "StringHelper.hpp"
#include "AtomicModelDataBlock.hpp"
#include "Logger.hpp"

#include <fstream>
#include <cstring>
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
    (void)model_object;
    (void)stream;
    (void)par;
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
