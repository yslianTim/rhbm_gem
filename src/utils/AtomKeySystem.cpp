#include <rhbm_gem/utils/AtomKeySystem.hpp>
#include <rhbm_gem/utils/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/StringHelper.hpp>
#include <rhbm_gem/utils/Logger.hpp>

#include <algorithm>

const AtomKey AtomKeySystem::k_dynamic_base{ 10000 };
const AtomKey AtomKeySystem::k_max_key{ std::numeric_limits<AtomKey>::max() };

AtomKeySystem::AtomKeySystem() :
    m_next_dynamic_key{ k_dynamic_base }
{
    const auto & build_in_spot_map{ ChemicalDataHelper::GetSpotMap() };
    for (const auto & [atom_id, spot] : build_in_spot_map)
    {
        auto atom_key{ static_cast<AtomKey>(spot) };
        m_id_to_key_map.emplace(atom_id, atom_key);
        m_key_to_id_map.emplace(atom_key, std::string{atom_id});
    }
}

AtomKeySystem::AtomKeySystem(const AtomKeySystem & other) :
    m_next_dynamic_key{ other.m_next_dynamic_key },
    m_id_to_key_map{ other.m_id_to_key_map },
    m_key_to_id_map{ other.m_key_to_id_map }
{
}

AtomKeySystem::~AtomKeySystem()
{
}

void AtomKeySystem::RegisterAtom(const std::string & atom_id)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(atom_id) != m_id_to_key_map.end()) return;
    if (m_next_dynamic_key >= k_max_key)
    {
        Logger::Log(LogLevel::Warning,
            "AtomKeySystem::RegisterComponent() - Atom key overflow for " + atom_id
            + ", avoiding register new atom key duw to maximum key reached."
        );
        return;
    }
    AtomKey new_atom_key{ m_next_dynamic_key++ };
    m_id_to_key_map[atom_id] = new_atom_key;
    m_key_to_id_map[new_atom_key] = atom_id;
}

void AtomKeySystem::RegisterAtom(const std::string & atom_id, AtomKey atom_key)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(atom_id) != m_id_to_key_map.end()) return;
    m_id_to_key_map[atom_id] = atom_key;
    m_key_to_id_map[atom_key] = atom_id;
}

AtomKey AtomKeySystem::GetAtomKey(const std::string & atom_id)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(atom_id) == m_id_to_key_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "AtomKeySystem::GetAtomKey() - Unknown atom id: " + atom_id);
        return static_cast<AtomKey>(0);
    }
    return m_id_to_key_map.at(atom_id);
}

std::string AtomKeySystem::GetAtomId(AtomKey atom_key)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_key_to_id_map.find(atom_key) == m_key_to_id_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "AtomKeySystem::GetAtomId() - Unknown atom key: "+ std::to_string(atom_key));
        return "UNK";
    }
    return m_key_to_id_map.at(atom_key);
}

bool AtomKeySystem::IsBuildInAtom(const std::string & atom_id) const
{
    if (m_id_to_key_map.find(atom_id) == m_id_to_key_map.end()) return false;
    return m_id_to_key_map.at(atom_id) < k_dynamic_base;
}

bool AtomKeySystem::IsBuildInAtom(AtomKey atom_key) const
{
    return atom_key < k_dynamic_base;
}
