#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <unordered_map>

enum class Residue : uint16_t;
enum class Element : uint16_t;
enum class Remoteness : uint8_t;
enum class Branch : uint8_t;

class AminoAcidInfoHelper
{
    static const std::unordered_map<Residue, size_t> m_atom_count_map;
    static const std::unordered_map<Residue, std::vector<Element>> m_element_map;
    static const std::unordered_map<Residue, std::vector<Remoteness>> m_remoteness_map;
    static const std::unordered_map<Residue, std::vector<Branch>> m_branch_map;
    static const std::unordered_map<Residue, std::vector<double>> m_buried_partial_charge_map;

public:
    AminoAcidInfoHelper(void) = default;
    ~AminoAcidInfoHelper() = default;

    static size_t GetAtomCount(Residue residue);
    static size_t GetAtomCount(int residue);
    static double GetPartialCharge(Residue residue, Element element, Remoteness remoteness, Branch branch, bool verbose=0);

};