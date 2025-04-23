#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <tuple>
#include <unordered_map>

enum class Residue : uint16_t;
enum class Element : uint16_t;
enum class Remoteness : uint8_t;
enum class Branch : uint8_t;
enum class Structure : uint8_t;

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

    static const std::unordered_map<Element, int> m_element_color_map;
    static const std::unordered_map<Residue, int> m_residue_color_map;
    static const std::unordered_map<Element, int> m_element_marker_map;
    static const std::unordered_map<Residue, int> m_residue_marker_map;

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

    static short GetDisplayColor(Element element);
    static short GetDisplayColor(Residue residue);
    static short GetDisplayMarker(Element element);
    static short GetDisplayMarker(Residue residue);
private:
 
};