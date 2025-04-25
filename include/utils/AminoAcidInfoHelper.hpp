#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>

enum class Residue : uint16_t;
enum class Element : uint16_t;
enum class Remoteness : uint8_t;
enum class Branch : uint8_t;
enum class Structure : uint8_t;

class AminoAcidInfoHelper
{
    static const std::unordered_map<Residue, std::vector<Element>> m_element_map;
    static const std::unordered_map<Residue, std::vector<Remoteness>> m_remoteness_map;
    static const std::unordered_map<Residue, std::vector<Branch>> m_branch_map;
    static const std::unordered_map<Residue, std::vector<double>> m_buried_partial_charge_map;
    static const std::unordered_map<Residue, std::vector<double>> m_helix_partial_charge_map;

public:
    AminoAcidInfoHelper(void) = default;
    ~AminoAcidInfoHelper() = default;

    static size_t GetAtomCount(Residue residue);
    static size_t GetAtomCount(int residue);
    static double GetPartialCharge(
        Residue residue, Element element, Remoteness remoteness, Branch branch,
        Structure structure, bool verbose=0);
    static const std::vector<double> & GetPartialChargeList(Residue residue, Structure structure);

};