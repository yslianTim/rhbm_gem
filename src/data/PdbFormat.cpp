#include "PdbFormat.hpp"
#include "AtomObject.hpp"
#include "StringHelper.hpp"

#include <fstream>
#include <iostream>
#include <cstring>
#include <stdexcept>

#define _CRT_SECURE_NO_WARNINGS // To disable deprecation warnings for sscanf

PdbFormat::PdbFormat(void)
{

}

PdbFormat::~PdbFormat()
{

}

void PdbFormat::LoadHeader(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadHeader failed!");
    }
}

void PdbFormat::PrintHeader(void) const
{

}

void PdbFormat::LoadDataArray(const std::string & filename)
{
    LoadAtomSiteData(filename);
}

void PdbFormat::LoadAtomSiteData(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadAtomSiteData failed!");
    }

    std::string line;
    while (std::getline(infile, line))
    {
        char * buffer{ line.data() };
        std::string header_str{ line.substr(0, header_string_length) };
        if (header_str == "ENDMDL") break; // end of reading file for multi-model PDB
        switch (MapToHeaderType(header_str))
        {
        case PDB_HEADER::ATOM:
            ScanAtomEntry(buffer, false);
            break;
        case PDB_HEADER::HETATM:
            ScanAtomEntry(buffer, true);
            break;
        default:
            break;
        }
    }
}

void PdbFormat::ScanAtomEntry(char * line, bool is_special)
{
    std::shared_ptr<AtomSite> atom{ std::make_shared<AtomSite>() };
    std::sscanf(&line[static_cast<int>(ATOM::SERIAL_ID)],    "%d",  &atom->serial_id);
    std::sscanf(&line[static_cast<int>(ATOM::ATOM_NAME)],    "%4c", &atom->atom_name[0]);
    std::sscanf(&line[static_cast<int>(ATOM::INDICATOR)],    "%c",  &atom->indicator);
    std::sscanf(&line[static_cast<int>(ATOM::RESIDUE_NAME)], "%3c", &atom->residue_name[0]);
    std::sscanf(&line[static_cast<int>(ATOM::CHAIN_ID)],     "%c",  &atom->chain_id);
    std::sscanf(&line[static_cast<int>(ATOM::RESIDUE_ID)],   "%d",  &atom->residue_id);
    std::sscanf(&line[static_cast<int>(ATOM::CODE)],         "%c",  &atom->code);
    std::sscanf(&line[static_cast<int>(ATOM::POSITION_X)],   "%8f", &atom->position_x);
    std::sscanf(&line[static_cast<int>(ATOM::POSITION_Y)],   "%8f", &atom->position_y);
    std::sscanf(&line[static_cast<int>(ATOM::POSITION_Z)],   "%8f", &atom->position_z);
    std::sscanf(&line[static_cast<int>(ATOM::OCCUPANCY)],    "%8f", &atom->occupancy);
    std::sscanf(&line[static_cast<int>(ATOM::TEMPERATURE)],  "%8f", &atom->temperature);
    std::sscanf(&line[static_cast<int>(ATOM::SEGMENT_ID)],   "%4c", &atom->segment_id[0]);
    std::sscanf(&line[static_cast<int>(ATOM::ELEMENT)],      "%2c", &atom->element[0]);
    std::sscanf(&line[static_cast<int>(ATOM::CHARGE)],       "%2c", &atom->charge[0]);
    BuildAtomObject(atom, is_special);
}

void PdbFormat::BuildAtomObject(std::any atom_info, bool is_special_atom)
{
    try
    {
        auto atom{ std::any_cast<std::shared_ptr<AtomSite>>(atom_info) };
        auto atom_object{ std::make_unique<AtomObject>() };
        auto atom_name{ StringHelper::ConvertCharArrayToString(atom->atom_name) };
        auto element_name{ StringHelper::ConvertCharArrayToString(atom->element) };
        auto indicator{ (atom->indicator == ' ') ? "." : std::string(1, atom->indicator) };
        atom_object->SetElement(element_name);
        atom_object->SetRemoteness(StringHelper::ExtractCharAsString(atom_name, element_name.size()));
        atom_object->SetBranch(StringHelper::ExtractCharAsString(atom_name, element_name.size()+1));
        atom_object->SetResidue(StringHelper::ConvertCharArrayToString(atom->residue_name));
        atom_object->SetIndicator(indicator);
        atom_object->SetResidueID(atom->residue_id);
        atom_object->SetSerialID(atom->serial_id);
        atom_object->SetChainID(std::string(1, atom->chain_id));
        atom_object->SetPosition(atom->position_x, atom->position_y, atom->position_z);
        atom_object->SetOccupancy(atom->occupancy);
        atom_object->SetTemperature(atom->temperature);
        atom_object->SetSpecialAtomFlag(is_special_atom);
        //atom_object->SetCharge();
        m_atom_object_list.emplace_back(std::move(atom_object));
    }
    catch (const std::bad_any_cast &)
    {
        std::cerr << "Error: bad any cast in BuildAtomObject" << std::endl;
    }
    catch (const std::exception & e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

std::vector<std::unique_ptr<AtomObject>> PdbFormat::GetAtomObjectList(void)
{
    return std::move(m_atom_object_list);
}

std::string PdbFormat::GetPdbID(void) const
{
    return m_model_id;
}

std::string PdbFormat::GetEmdID(void) const
{
    return m_map_id;
}

PdbFormat::PDB_HEADER PdbFormat::MapToHeaderType(const std::string & name) const
{
    try
    {
        header_map.at(name);
    }
    catch(const std::exception & except)
    {
        std::cerr << except.what() << std::endl;
        return PDB_HEADER::UNK;
    }
    return header_map.at(name);
}