#include "BondKeySystem.hpp"
#include "GlobalEnumClass.hpp"
#include "ChemicalDataHelper.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

#include <algorithm>

const BondKey BondKeySystem::k_dynamic_base{ 10000 };
const BondKey BondKeySystem::k_max_key{ std::numeric_limits<BondKey>::max() };

BondKeySystem::BondKeySystem(void) :
    m_next_dynamic_key{ k_dynamic_base }
{
    Logger::Log(LogLevel::Debug, "BondKeySystem::BondKeySystem() called");
    for (const auto & [bond_id, link] : ChemicalDataHelper::GetLinkMap())
    {
        auto bond_key{ static_cast<BondKey>(link) };
        m_id_to_key_map.emplace(bond_id, bond_key);
        m_key_to_id_map.emplace(bond_key, std::string{bond_id});
        m_veto_bond_id_set.emplace(BuildReverseBondIdFromBondId(std::string{bond_id}));
    }
}

BondKeySystem::BondKeySystem(const BondKeySystem & other) :
    m_next_dynamic_key{ other.m_next_dynamic_key },
    m_id_to_key_map{ other.m_id_to_key_map },
    m_key_to_id_map{ other.m_key_to_id_map }
{
    Logger::Log(LogLevel::Debug, "BondKeySystem::BondKeySystem() called");
}

BondKeySystem::~BondKeySystem()
{
    Logger::Log(LogLevel::Debug, "BondKeySystem::~BondKeySystem() called");
}

void BondKeySystem::RegisterBond(
    const std::string & bond_id)
{
    Logger::Log(LogLevel::Debug, "BondKeySystem::RegisterBond() called");
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(bond_id) != m_id_to_key_map.end()) return;
    if (m_next_dynamic_key >= k_max_key)
    {
        Logger::Log(LogLevel::Warning,
            "BondKeySystem::RegisterComponent() - Bond key overflow for " + bond_id
            + ", avoiding register new atom key duw to maximum key reached."
        );
        return;
    }
    BondKey new_bond_key{ m_next_dynamic_key++ };
    m_id_to_key_map[bond_id] = new_bond_key;
    m_key_to_id_map[new_bond_key] = bond_id;
}

void BondKeySystem::RegisterBond(
    const std::string & atom_id_1, const std::string & atom_id_2)
{
    auto bond_id{ BuildBondIdFromAtomIdPair(atom_id_1, atom_id_2) };
    RegisterBond(bond_id);
}

void BondKeySystem::RegisterBond(
    const std::string & bond_id, BondKey bond_key)
{
    Logger::Log(LogLevel::Debug, "BondKeySystem::RegisterBond() called");
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(bond_id) != m_id_to_key_map.end()) return;
    m_id_to_key_map[bond_id] = bond_key;
    m_key_to_id_map[bond_key] = bond_id;
}

void BondKeySystem::RegisterBond(
    const std::string & atom_id_1, const std::string & atom_id_2, BondKey bond_key)
{
    auto bond_id{ BuildBondIdFromAtomIdPair(atom_id_1, atom_id_2) };
    RegisterBond(bond_id, bond_key);
}

BondKey BondKeySystem::GetBondKey(
    const std::string & bond_id)
{
    static std::unordered_map<std::string, int> unknown_bond_id_count_map;
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(bond_id) == m_id_to_key_map.end())
    {
        if (unknown_bond_id_count_map.find(bond_id) == unknown_bond_id_count_map.end())
        {
            Logger::Log(LogLevel::Warning,
                "BondKeySystem::GetBondKey() - Unknown bond id: [" + bond_id + "]");
            unknown_bond_id_count_map[bond_id] = 1;
        }
        else
        {
            unknown_bond_id_count_map[bond_id]++;
        }
        return static_cast<BondKey>(0);
    }
    return m_id_to_key_map.at(bond_id);
}

BondKey BondKeySystem::GetBondKey(
    const std::string & atom_id_1, const std::string & atom_id_2)
{
    static std::unordered_map<std::string, int> unknown_bond_id_count_map;
    auto bond_id{ BuildBondIdFromAtomIdPair(atom_id_1, atom_id_2) };
    return GetBondKey(bond_id);
}

std::string BondKeySystem::GetBondId(BondKey bond_key)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_key_to_id_map.find(bond_key) == m_key_to_id_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "BondKeySystem::GetBondId() - Unknown bond key: "+ std::to_string(bond_key));
        return "UNK";
    }
    return m_key_to_id_map.at(bond_key);
}

bool BondKeySystem::IsRegistedBond(
    const std::string & atom_id_1, const std::string & atom_id_2) const
{
    auto bond_id{ BuildBondIdFromAtomIdPair(atom_id_1, atom_id_2) };
    if (m_id_to_key_map.find(bond_id) == m_id_to_key_map.end()) return false;
    return true;
}

bool BondKeySystem::IsBuildInBond(
    const std::string & atom_id_1, const std::string & atom_id_2) const
{
    auto bond_id{ BuildBondIdFromAtomIdPair(atom_id_1, atom_id_2) };
    if (m_id_to_key_map.find(bond_id) == m_id_to_key_map.end()) return false;
    return m_id_to_key_map.at(bond_id) < k_dynamic_base;
}

bool BondKeySystem::IsBuildInBond(BondKey bond_key) const
{
    return bond_key < k_dynamic_base;
}

bool BondKeySystem::IsReverseBond(
    const std::string & atom_id_1, const std::string & atom_id_2) const
{
    auto bond_id{ BuildBondIdFromAtomIdPair(atom_id_1, atom_id_2) };
    return m_veto_bond_id_set.find(bond_id) != m_veto_bond_id_set.end();
}

std::string BondKeySystem::BuildBondIdFromAtomIdPair(
    const std::string & atom_id_1, const std::string & atom_id_2) const
{
    return atom_id_1 +"_"+ atom_id_2;
}

std::string BondKeySystem::BuildReverseBondIdFromBondId(
    const std::string & bond_id) const
{
    auto atom_id_pair{ BuildAtomIdPairFromBondId(bond_id) };
    if (atom_id_pair.first == atom_id_pair.second) return ".";
    return atom_id_pair.second +"_"+ atom_id_pair.first;
}

std::pair<std::string, std::string> BondKeySystem::BuildAtomIdPairFromBondId(
    const std::string & bond_id) const
{
    auto pos{ bond_id.find("_") };
    if (pos == std::string::npos)
    {
        Logger::Log(LogLevel::Warning, 
            "BondKeySystem::BuildAtomIdPairFromBondId() - Invalid bond id: " + bond_id);
        return {"", ""};
    }
    auto atom_id_1{ bond_id.substr(0, pos) };
    auto atom_id_2{ bond_id.substr(pos+1) };
    return { atom_id_1, atom_id_2 };
}
