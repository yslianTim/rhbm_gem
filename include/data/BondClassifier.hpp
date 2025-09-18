#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>

#include "GlobalEnumClass.hpp"

class BondObject;
class BondClassifier
{
    static const size_t m_main_chain_member_count{ 4 };
    static const std::vector<short> m_main_chain_member_color_list;
    static const std::vector<short> m_main_chain_member_solid_marker_list;
    static const std::vector<short> m_main_chain_member_open_marker_list;
    static const std::vector<Bond> m_main_chain_member_bond_list;
    static const std::vector<std::string> m_main_chain_member_label_list;
    static const std::vector<std::string> m_main_chain_member_title_list;

public:
    BondClassifier(void);
    ~BondClassifier();

    static bool IsMainChainMember(Bond bond, size_t & main_chain_member_id);
    static size_t GetMainChainMemberCount(void);
    static Bond GetMainChainBond(size_t id);
    static short GetMainChainMemberColor(size_t id);
    static short GetMainChainMemberSolidMarker(size_t id);
    static short GetMainChainMemberOpenMarker(size_t id);
    static const std::string & GetMainChainMemberLabel(size_t id);
    static const std::string & GetMainChainMemberTitle(size_t id);

    static GroupKey GetGroupKeyInClass(const BondObject * bond_object, const std::string & class_key);

    GroupKey GetMainChainSimpleBondClassGroupKey(size_t id) const;
    GroupKey GetMainChainComponentBondClassGroupKey(size_t id, Residue residue) const;
    std::vector<GroupKey> GetMainChainComponentBondClassGroupKeyList(size_t id) const;

private:
    static bool IsValidMainChainMemberID(size_t id);

};
