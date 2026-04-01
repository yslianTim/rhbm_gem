#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>

namespace rhbm_gem {

class AtomObject;
class AtomClassifier
{
    static const size_t m_nucleotide_main_chain_member_start_id{ 4 };
    static const size_t m_nucleotide_main_chain_member_count{ 3 + 9 };
    static const size_t m_main_chain_member_count{ 4 + m_nucleotide_main_chain_member_count };
    static const std::vector<short> m_main_chain_member_color_list;
    static const std::vector<short> m_main_chain_member_solid_marker_list;
    static const std::vector<short> m_main_chain_member_open_marker_list;
    static const std::unordered_map<Spot, short> m_main_chain_color_map;
    static const std::unordered_map<Spot, short> m_main_chain_solid_marker_map;
    static const std::unordered_map<Spot, short> m_main_chain_open_marker_map;
    static const std::unordered_map<Spot, short> m_main_chain_line_style_map;
    static const std::unordered_map<Spot, std::string> m_main_chain_spot_label_map;
    static const std::unordered_map<Structure, short> m_main_chain_structure_color_map;
    static const std::unordered_map<Structure, short> m_main_chain_structure_marker_map;
    static const std::unordered_map<Structure, short> m_main_chain_structure_line_style_map;
    static const std::unordered_map<Structure, std::string> m_main_chain_structure_label_map;
    static const std::vector<Element> m_main_chain_member_element_list;
    static const std::vector<Spot> m_main_chain_member_spot_list;
    static const std::vector<std::string> m_main_chain_member_label_list;
    static const std::vector<std::string> m_main_chain_member_title_list;

public:
    static bool IsMainChainMember(Spot spot, size_t & main_chain_member_id);
    static size_t GetMainChainMemberCount();
    static Element GetMainChainElement(size_t member_id);
    static Spot GetMainChainSpot(size_t member_id);
    static Spot GetNucleotideMainChainSpot(size_t member_id);
    static short GetMainChainElementColor(size_t member_id);
    static short GetMainChainSpotColor(Spot spot);
    static short GetMainChainStructureColor(Structure structure);
    static short GetNucleotideMainChainElementColor(size_t member_id);
    static short GetMainChainElementSolidMarker(size_t member_id);
    static short GetMainChainSpotSolidMarker(Spot spot);
    static short GetMainChainStructureMarker(Structure structure);
    static short GetMainChainElementOpenMarker(size_t member_id);
    static short GetMainChainSpotOpenMarker(Spot spot);
    static short GetMainChainSpotLineStyle(Spot spot);
    static short GetMainChainStructureLineStyle(Structure structure);
    static const std::string & GetMainChainElementLabel(size_t member_id);
    static const std::string & GetMainChainSpotLabel(Spot spot);
    static const std::string & GetMainChainStructureLabel(Structure structure);
    static const std::string & GetMainChainElementTitle(size_t member_id);

    static GroupKey GetGroupKeyInClass(const AtomObject * atom_object, const std::string & class_key);
    static GroupKey GetGroupKeyInClass(ComponentKey component_key, AtomKey atom_key);
    static GroupKey GetGroupKeyInClass(Structure structure, ComponentKey component_key, AtomKey atom_key);

    static GroupKey GetMainChainSimpleAtomClassGroupKey(size_t member_id);
    static GroupKey GetMainChainComponentAtomClassGroupKey(size_t member_id, Residue residue);
    static GroupKey GetMainChainStructureAtomClassGroupKey(size_t member_id, Structure structure, Residue residue);
    static std::vector<GroupKey> GetMainChainComponentAtomClassGroupKeyList(size_t member_id);
    static std::vector<GroupKey> GetMainChainStructureAtomClassGroupKeyList(size_t member_id, Structure structure);

    static GroupKey GetNucleotideMainChainComponentAtomClassGroupKey(size_t member_id, Residue residue);
    static std::vector<GroupKey> GetNucleotideMainChainComponentAtomClassGroupKeyList(Residue component_id);

private:
    static bool IsValidMainChainMemberID(size_t member_id);

};

} // namespace rhbm_gem
