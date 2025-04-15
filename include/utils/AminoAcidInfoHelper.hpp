#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>
#include "AtomicInfoHelper.hpp"

class AminoAcidInfoHelper
{

    inline static std::unordered_map<Residue, int> amino_acid_atom_count_map
    {
        {Residue::ALA,  5}, {Residue::ARG, 10}, {Residue::ASN,  7}, {Residue::ASP,  7},
        {Residue::CYS,  6}, {Residue::GLN,  8}, {Residue::GLU,  8}, {Residue::GLY,  4},
        {Residue::HIS,  8}, {Residue::ILE,  7}, {Residue::LEU,  7}, {Residue::LYS,  9},
        {Residue::MET,  8}, {Residue::PHE,  9}, {Residue::PRO,  7}, {Residue::SER,  6},
        {Residue::THR,  6}, {Residue::TRP, 10}, {Residue::TYR, 10}, {Residue::VAL,  6}
    };

public:
    AminoAcidInfoHelper(void) = default;
    ~AminoAcidInfoHelper() = default;

};