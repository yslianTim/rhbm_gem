#pragma once

#include <iostream>
#include <string>
#include <string_view>
#include <sstream>
#include <tuple>
#include <exception>
#include <unordered_map>

enum class Residue : int
{
    ALA =  0, ARG =  1, ASN =  2, ASP =  3, CYS =  4,
    GLN =  5, GLU =  6, GLY =  7, HIS =  8, ILE =  9,
    LEU = 10, LYS = 11, MET = 12, PHE = 13, PRO = 14,
    SER = 15, THR = 16, TRP = 17, TYR = 18, VAL = 19,
    CSX = 20,
    HOH = 99,
    NA = 911, MG = 912, CL = 917, CA = 920, FE = 926, FE2 = 826, ZN = 930,
    UNK = -1
};

enum class Element : int
{
    HYDROGEN = 1, CARBON = 6, NITROGEN = 7, OXYGEN = 8,
    SODIUM = 11, MAGNESIUM = 12, SULFUR = 16, CHLORINE = 17,
    CALCIUM = 20, IRON = 26, ZINC = 30,
    UNK = -1
};

enum class Remoteness : int
{
    NONE = 0, ALPHA = 1, BETA = 2, GAMMA = 3, DELTA = 4, EPSILON = 5, ZETA = 6, ETA = 7,
    ONE = 11, TWO = 12, THREE = 13, FOUR = 14, FIVE = 15,
    EXTRA = 99,
    UNK = -1
};

enum class Branch : int
{
    NONE = 0, ONE = 1, TWO = 2, THREE = 3,
    TERMINAL = 99,
    UNK = -1
};

class AtomicInfoHelper
{
    inline static std::unordered_map<std::string_view, Residue> residue_map
    {
        {"ALA", Residue::ALA}, {"ARG", Residue::ARG}, {"ASN", Residue::ASN}, {"ASP", Residue::ASP},
        {"CYS", Residue::CYS}, {"GLN", Residue::GLN}, {"GLU", Residue::GLU}, {"GLY", Residue::GLY},
        {"HIS", Residue::HIS}, {"ILE", Residue::ILE}, {"LEU", Residue::LEU}, {"LYS", Residue::LYS},
        {"MET", Residue::MET}, {"PHE", Residue::PHE}, {"PRO", Residue::PRO}, {"SER", Residue::SER},
        {"THR", Residue::THR}, {"TRP", Residue::TRP}, {"TYR", Residue::TYR}, {"VAL", Residue::VAL},
        {"HOH", Residue::HOH}, {"NA",  Residue::NA},  {"MG",  Residue::MG},  {"CA",  Residue::CA},
        {"ZN",  Residue::ZN},  {"FE",  Residue::FE},  {"FE2", Residue::FE2}, {"CSX", Residue::CSX},
        {"CL",  Residue::CL},  {"UNK", Residue::UNK}
    };

    inline static std::unordered_map<std::string_view, Element> element_map
    {
        {"H",  Element::HYDROGEN},
        {"C",  Element::CARBON},    {"N",  Element::NITROGEN},
        {"O",  Element::OXYGEN},    {"S",  Element::SULFUR},
        {"CA", Element::CALCIUM},   {"ZN", Element::ZINC},    {"NA", Element::SODIUM},
        {"MG", Element::MAGNESIUM}, {"FE", Element::IRON},    {"CL", Element::CHLORINE},
        {"UNK",Element::UNK}
    };

    inline static std::unordered_map<std::string_view, Remoteness> remoteness_map
    {
        {" ", Remoteness::NONE}, {"A", Remoteness::ALPHA}, {"B", Remoteness::BETA},
        {"G", Remoteness::GAMMA},{"D", Remoteness::DELTA}, {"E", Remoteness::EPSILON},
        {"Z", Remoteness::ZETA}, {"H", Remoteness::ETA},
        {"1", Remoteness::ONE},  {"2", Remoteness::TWO},   {"3", Remoteness::THREE},
        {"4", Remoteness::FOUR}, {"5", Remoteness::FIVE},
        {"X", Remoteness::EXTRA},{"UNK",Remoteness::UNK}
    };

    inline static std::unordered_map<std::string_view, Branch> branch_map
    {
        {" ", Branch::NONE}, {"1", Branch::ONE}, {"2", Branch::TWO}, {"3", Branch::THREE},
        {"T", Branch::TERMINAL}, {"UNK", Branch::UNK}
    };

public:
    AtomicInfoHelper(void) = default;
    ~AtomicInfoHelper() = default;
    struct ResidueTag {};
    struct ElementTag {};
    struct RemotenessTag {};
    struct BranchTag {};

    template <typename Tag>
    static int AtomInfoMapping(const std::string & name)
    {
        const auto & mapping{ GetMapping(Tag{}) };
        CheckEnumTypeName(name, mapping, typeid(Tag).name());
        return static_cast<int>(mapping.at(name));
    }

    static const std::unordered_map<std::string_view, Residue> & GetMapping(const ResidueTag &)
    {
        return residue_map;
    }

    static const std::unordered_map<std::string_view, Element> & GetMapping(const ElementTag &)
    {
        return element_map;
    }

    static const std::unordered_map<std::string_view, Remoteness> & GetMapping(const RemotenessTag &)
    {
        return remoteness_map;
    }

    static const std::unordered_map<std::string_view, Branch> & GetMapping(const BranchTag &)
    {
        return branch_map;
    }

private:
    template <typename EnumType>
    static void CheckEnumTypeName(const std::string & name,
                                  const std::unordered_map<std::string_view, EnumType> & mapping,
                                  const std::string & mapping_type)
    {
        if (mapping.find(name) == mapping.end())
        {
            std::ostringstream error_message;
            error_message << "Invalid " << mapping_type << " name: " << name;
            throw std::runtime_error(error_message.str());
        }
    }
};