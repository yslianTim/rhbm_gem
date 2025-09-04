#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include <unordered_map>

#include "GlobalEnumClass.hpp"

class AminoAcidInfoHelper
{
    static const std::unordered_map<Residue, std::vector<Spot>> m_spot_map;
    static const std::unordered_map<Residue, std::vector<Element>> m_element_map;
    static const std::unordered_map<Residue, std::vector<double>> m_amber95_partial_charge_map;
    static const std::unordered_map<Residue, std::vector<double>> m_buried_partial_charge_map;
    static const std::unordered_map<Residue, std::vector<double>> m_helix_partial_charge_map;
    static const std::unordered_map<Residue, std::vector<double>> m_sheet_partial_charge_map;

public:
    AminoAcidInfoHelper(void) = default;
    ~AminoAcidInfoHelper() = default;

    static size_t GetAtomCount(Residue residue);
    static size_t GetAtomCount(int residue);
    static double GetPartialCharge(
        Residue residue, Spot spot, Structure structure,
        bool use_amber_table=false, bool verbose=false
    );
    static const std::vector<double> & GetPartialChargeList(Residue residue, Structure structure);
    static const std::vector<double> & GetPartialChargeListAmber(Residue residue);
    static const std::vector<Spot> & GetSpotList(Residue residue);

};
