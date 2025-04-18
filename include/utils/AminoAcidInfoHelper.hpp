#pragma once

#include <cstddef>
#include <vector>
#include <unordered_map>

enum class Element : int;
enum class Residue : int;
enum class Remoteness : int;
enum class Branch : int;

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