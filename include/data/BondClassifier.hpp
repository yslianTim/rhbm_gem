#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>

#include "GlobalEnumClass.hpp"
#include "ComponentKeySystem.hpp"
#include "BondKeySystem.hpp"

class BondObject;
class BondClassifier
{
    static const size_t m_main_chain_member_count{ 4 };
    static const std::vector<short> m_main_chain_member_color_list;
    static const std::vector<short> m_main_chain_member_solid_marker_list;
    static const std::vector<short> m_main_chain_member_open_marker_list;
    static const std::vector<Link> m_main_chain_member_link_list;
    static const std::vector<std::string> m_main_chain_member_label_list;
    static const std::vector<std::string> m_main_chain_member_title_list;

public:
    BondClassifier(void);
    ~BondClassifier();

    static bool IsMainChainMember(Link link, size_t & main_chain_member_id);
    static size_t GetMainChainMemberCount(void);
    static Link GetMainChainLink(size_t id);
    static short GetMainChainMemberColor(size_t id);
    static short GetMainChainMemberSolidMarker(size_t id);
    static short GetMainChainMemberOpenMarker(size_t id);
    static const std::string & GetMainChainMemberLabel(size_t id);
    static const std::string & GetMainChainMemberTitle(size_t id);

    static GroupKey GetGroupKeyInClass(const BondObject * bond_object, const std::string & class_key);
    static GroupKey GetGroupKeyInClass(ComponentKey component_key, BondKey bond_key);

    static GroupKey GetMainChainSimpleBondClassGroupKey(size_t id);
    GroupKey GetMainChainComponentBondClassGroupKey(size_t id, Residue residue) const;
    std::vector<GroupKey> GetMainChainComponentBondClassGroupKeyList(size_t id) const;

private:
    static bool IsValidMainChainMemberID(size_t id);

};
