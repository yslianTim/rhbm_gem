#pragma once

#include <iostream>
#include <cuchar>
#include <vector>
#include <unordered_map>
#include "AtomicInfoHelper.hpp"

class AminoAcidInfoHelper
{
    inline static std::unordered_map<Residue, size_t> amino_acid_atom_count_map
    {
        {Residue::ALA,  5}, {Residue::ARG, 10}, {Residue::ASN,  7}, {Residue::ASP,  7},
        {Residue::CYS,  6}, {Residue::GLN,  8}, {Residue::GLU,  8}, {Residue::GLY,  4},
        {Residue::HIS,  8}, {Residue::ILE,  7}, {Residue::LEU,  7}, {Residue::LYS,  9},
        {Residue::MET,  8}, {Residue::PHE,  9}, {Residue::PRO,  7}, {Residue::SER,  6},
        {Residue::THR,  6}, {Residue::TRP, 10}, {Residue::TYR, 10}, {Residue::VAL,  6}
    };

    inline static std::unordered_map<Residue, std::vector<Element>> amino_acid_element_map
    {
        {Residue::ALA, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::ARG, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::ASN, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::ASP, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::CYS, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::SULFUR}},
        {Residue::GLN, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::GLU, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::GLY, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN}},
        {Residue::HIS, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::ILE, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::LEU, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::LYS, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::MET, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::SULFUR}},
        {Residue::PHE, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::PRO, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::SER, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::THR, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::TRP, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::TYR, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}},
        {Residue::VAL, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                        Element::CARBON}}
    };

    inline static std::unordered_map<Residue, std::vector<Remoteness>> amino_acid_remoteness_map
    {
        {Residue::ALA, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::ARG, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::ASN, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::ASP, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::CYS, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::GLN, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::GLU, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::GLY, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE}},
        {Residue::HIS, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::ILE, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::LEU, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::LYS, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::MET, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::PHE, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::PRO, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::SER, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::THR, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::TRP, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::TYR, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}},
        {Residue::VAL, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                        Remoteness::BETA}}
    };

    inline static std::unordered_map<Residue, std::vector<Branch>> amino_acid_branch_map
    {
        {Residue::ALA, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::ARG, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::ASN, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::ASP, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::CYS, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::GLN, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::GLU, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::GLY, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE}},
        {Residue::HIS, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::ILE, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::LEU, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::LYS, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::MET, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::PHE, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::PRO, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::SER, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::THR, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::TRP, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::TYR, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}},
        {Residue::VAL, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                        Branch::NONE}}
    };

public:
    AminoAcidInfoHelper(void) = default;
    ~AminoAcidInfoHelper() = default;

    static size_t GetAtomCount(Residue residue) { return amino_acid_atom_count_map.at(residue); }
    static size_t GetAtomCount(int residue) { return amino_acid_atom_count_map.at(static_cast<Residue>(residue)); }

};