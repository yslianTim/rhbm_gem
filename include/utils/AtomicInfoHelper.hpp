#pragma once

#include <iostream>
#include <string>
#include <string_view>
#include <sstream>
#include <vector>
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

struct ResidueGroupClassifierTag {};
struct ElementGroupClassifierTag {};

template <typename Tag>
struct GroupKeyMapping;

template <>
struct GroupKeyMapping<ResidueGroupClassifierTag>
{
    using type = std::tuple<int, int, int, int, bool>;
};

template <>
struct GroupKeyMapping<ElementGroupClassifierTag>
{
    using type = std::tuple<int, int, bool>;
};

class AtomicInfoHelper
{
    inline static int m_standard_residue_count{ 20 };
    inline static std::string m_element_class_key{ "element_class" };
    inline static std::string m_residue_class_key{ "residue_class" };
    static const std::vector<int> m_standard_residue_list;
    static const std::vector<int> m_standard_element_list;
    static const std::vector<int> m_standard_remoteness_list;
    static const std::vector<int> m_standard_branch_list;

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
    static int GetAtomicNumber(Element element) { return m_atomic_number_map.at(element); }
    static int GetStandardResidueCount(void) { return m_standard_residue_count; }
    static const std::string & GetElementClassKey(void) { return m_element_class_key; }
    static const std::string & GetResidueClassKey(void) { return m_residue_class_key; }
    static const std::vector<int> & GetStandardResidueList(void) { return m_standard_residue_list; }
    static const std::vector<int> & GetStandardElementList(void) { return m_standard_element_list; }
    static const std::vector<int> & GetStandardRemotenessList(void) { return m_standard_remoteness_list; }
    static const std::vector<int> & GetStandardBranchList(void) { return m_standard_branch_list; }

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

    template <typename Tag>
    static std::string AtomLabelMapping(int value)
    {
        const auto & mapping{ GetLabelMapping(Tag{}) };
        using KeyType = typename std::remove_reference<decltype(mapping)>::type::key_type;
        KeyType key{ static_cast<KeyType>(value) };
        if (mapping.find(key) == mapping.end())
        {
            std::ostringstream error_message;
            error_message << "Invalid value: " << value << " for mapping type " << typeid(Tag).name();
            throw std::runtime_error(error_message.str());
        }
        return mapping.at(key);
    }

    static const std::unordered_map<Residue, std::string> & GetLabelMapping(const ResidueTag &)
    {
        return residue_label_map;
    }

    static const std::unordered_map<Element, std::string> & GetLabelMapping(const ElementTag &)
    {
        return element_label_map;
    }

    static const std::unordered_map<Remoteness, std::string> & GetLabelMapping(const RemotenessTag &)
    {
        return remoteness_label_map;
    }

    static const std::unordered_map<Branch, std::string> & GetLabelMapping(const BranchTag &)
    {
        return branch_label_map;
    }

    template <typename Tag>
    static int AtomColorMapping(int value)
    {
        const auto & mapping{ GetColorMapping(Tag{}) };
        using KeyType = typename std::remove_reference<decltype(mapping)>::type::key_type;
        KeyType key{ static_cast<KeyType>(value) };
        if (mapping.find(key) == mapping.end())
        {
            std::ostringstream error_message;
            error_message << "Invalid value: " << value << " for mapping type " << typeid(Tag).name();
            throw std::runtime_error(error_message.str());
        }
        return mapping.at(key);
    }

    static const std::unordered_map<Element, int> & GetColorMapping(const ElementTag &)
    {
        return element_color_map;
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