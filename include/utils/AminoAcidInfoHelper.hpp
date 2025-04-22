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
    static const std::unordered_map<Residue, size_t> atom_count_map;
    static const std::unordered_map<Residue, std::vector<Element>> element_map;
    static const std::unordered_map<Residue, std::vector<Remoteness>> remoteness_map;
    static const std::unordered_map<Residue, std::vector<Branch>> branch_map;
    static const std::unordered_map<Residue, std::vector<double>> buried_partial_charge_map;

public:
    AminoAcidInfoHelper(void) = default;
    ~AminoAcidInfoHelper() = default;

    static size_t GetAtomCount(Residue residue) { return atom_count_map.at(residue); }
    static size_t GetAtomCount(int residue) { return atom_count_map.at(static_cast<Residue>(residue)); }
    static double GetPartialCharge(Residue residue, Element element, Remoteness remoteness, Branch branch);

};