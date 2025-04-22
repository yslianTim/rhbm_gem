#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <tuple>
#include <unordered_map>

#ifndef UINT16_MAX
#define UINT16_MAX 65535
#endif

#ifndef UINT8_MAX
#define UINT8_MAX 255
#endif

enum class Residue : uint16_t
{
    ALA =  0, ARG =  1, ASN =  2, ASP =  3, CYS =  4,
    GLN =  5, GLU =  6, GLY =  7, HIS =  8, ILE =  9,
    LEU = 10, LYS = 11, MET = 12, PHE = 13, PRO = 14,
    SER = 15, THR = 16, TRP = 17, TYR = 18, VAL = 19,
    CSX = 20,
    HOH = 99,
    NA = 911, MG = 912, CL = 917, CA = 920, FE = 926, FE2 = 826, ZN = 930,
    UNK = UINT16_MAX
};

enum class Element : uint16_t
{
    HYDROGEN = 1, CARBON = 6, NITROGEN = 7, OXYGEN = 8,
    SODIUM = 11, MAGNESIUM = 12, SULFUR = 16, CHLORINE = 17,
    CALCIUM = 20, IRON = 26, ZINC = 30,
    UNK = UINT16_MAX
};

enum class Remoteness : uint8_t
{
    NONE = 0, ALPHA = 1, BETA = 2, GAMMA = 3, DELTA = 4, EPSILON = 5, ZETA = 6, ETA = 7,
    ONE = 11, TWO = 12, THREE = 13, FOUR = 14, FIVE = 15,
    EXTRA = 99,
    UNK = UINT8_MAX
};

enum class Branch : uint8_t
{
    NONE = 0, ONE = 1, TWO = 2, THREE = 3,
    TERMINAL = 254,
    UNK = UINT8_MAX
};

enum class Structure : uint8_t
{
    FREE = 0, HELIX = 1, SHEET = 2,
    UNK = UINT8_MAX
};

struct ResidueGroupClassifierTag {};
struct ElementGroupClassifierTag {};

template <typename Tag>
struct GroupKeyMapping;

template <>
struct GroupKeyMapping<ResidueGroupClassifierTag>
{
    using type = std::tuple<Residue, Element, Remoteness, Branch, bool>;
};

template <>
struct GroupKeyMapping<ElementGroupClassifierTag>
{
    using type = std::tuple<Element, Remoteness, bool>;
};

class AtomicInfoHelper
{
    inline static int m_standard_residue_count{ 20 };
    inline static std::string m_element_class_key{ "element_class" };
    inline static std::string m_residue_class_key{ "residue_class" };
    static const std::vector<Residue> m_standard_residue_list;
    static const std::vector<Element> m_standard_element_list;
    static const std::vector<Remoteness> m_standard_remoteness_list;
    static const std::vector<Branch> m_standard_branch_list;

    static const std::unordered_map<Element, int> m_atomic_number_map;
    static const std::unordered_map<std::string_view, Residue> residue_map;
    static const std::unordered_map<std::string_view, Element> element_map;
    static const std::unordered_map<std::string_view, Remoteness> remoteness_map;
    static const std::unordered_map<std::string_view, Branch> branch_map;
    static const std::unordered_map<Residue, std::string> residue_label_map;
    static const std::unordered_map<Element, std::string> element_label_map;
    static const std::unordered_map<Remoteness, std::string> remoteness_label_map;
    static const std::unordered_map<Branch, std::string> branch_label_map;

    static const std::unordered_map<Element, int> element_color_map;


public:
    AtomicInfoHelper(void) = default;
    ~AtomicInfoHelper() = default;
    static int GetAtomicNumber(Element element);
    static int GetStandardResidueCount(void);
    static const std::string & GetElementClassKey(void);
    static const std::string & GetResidueClassKey(void);
    static const std::vector<Residue> & GetStandardResidueList(void);
    static const std::vector<Element> & GetStandardElementList(void);
    static const std::vector<Remoteness> & GetStandardRemotenessList(void);
    static const std::vector<Branch> & GetStandardBranchList(void);
    static bool IsStandardResidue(Residue residue);

    static Residue GetResidueFromString(const std::string & name);
    static Element GetElementFromString(const std::string & name);
    static Remoteness GetRemotenessFromString(const std::string & name);
    static Branch GetBranchFromString(const std::string & name);

    static const std::string & GetLabel(Residue residue);
    static const std::string & GetLabel(Element element);
    static const std::string & GetLabel(Remoteness remoteness);
    static const std::string & GetLabel(Branch branch);

    static int GetDisplayColor(Element element);
private:
 
};