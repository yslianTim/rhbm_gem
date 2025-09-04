#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>

#include "GlobalEnumClass.hpp"

class AtomicInfoHelper
{
    static const std::vector<std::string> m_group_class_key_list;
    static const std::vector<Residue> m_standard_residue_list;
    static const std::vector<Element> m_standard_element_list;

    static const std::unordered_map<int, Element> m_atomic_number_to_element_map;
    static const std::unordered_map<std::string_view, Residue> m_residue_map;
    static const std::unordered_map<std::string_view, Element> m_element_map;
    static const std::unordered_map<std::string_view, Spot> m_spot_map;
    static const std::unordered_map<std::string_view, Structure> m_structure_map;
    static const std::unordered_map<std::string_view, Entity> m_entity_map;

    static const std::unordered_map<Residue, std::string> m_residue_label_map;
    static const std::unordered_map<Element, std::string> m_element_label_map;

    static const std::unordered_map<Element, int> m_element_color_map;
    static const std::unordered_map<Residue, int> m_residue_color_map;
    static const std::unordered_map<Element, int> m_element_marker_map;
    static const std::unordered_map<Residue, int> m_residue_marker_map;

public:
    AtomicInfoHelper(void) = default;
    ~AtomicInfoHelper() = default;
    static int GetAtomicNumber(Element element);
    static size_t GetGroupClassCount(void);
    static size_t GetStandardResidueCount(void);
    static size_t GetElementCount(void);
    static const std::string & GetGroupClassKey(size_t class_id);
    static const std::string & GetElementClassKey(void);
    static const std::string & GetResidueClassKey(void);
    static const std::string & GetStructureClassKey(void);
    static const std::vector<Residue> & GetStandardResidueList(void);
    static const std::vector<Element> & GetStandardElementList(void);
    static const std::unordered_map<std::string_view, Residue> & GetResidueMap(void);
    static const std::unordered_map<std::string_view, Element> & GetElementMap(void);
    static const std::unordered_map<std::string_view, Spot> & GetSpotMap(void);
    static const std::unordered_map<Element, std::string> & GetElementLabelMap(void);
    static bool IsStandardElement(Element element);
    static bool IsStandardResidue(Residue residue);

    static Residue GetResidueFromString(const std::string & name, bool verbose=true);
    static Element GetElementFromString(const std::string & name);
    static Element GetElementFromAtomicNumber(int atomic_number);
    static Spot GetSpotFromString(const std::string & name, bool verbose=true);
    static Structure GetStructureFromString(const std::string & name);
    static Entity GetEntityFromString(const std::string & name);

    static const std::string & GetLabel(Residue residue);
    static const std::string & GetLabel(Element element);

    static short GetDisplayColor(Element element);
    static short GetDisplayColor(Residue residue);
    static short GetDisplayMarker(Element element);
    static short GetDisplayMarker(Residue residue);

};
