#include "AtomKeySystem.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomicInfoHelper.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

#include <algorithm>

const AtomKey AtomKeySystem::kDynamicBase{ static_cast<uint32_t>(Branch::TERMINAL) << 16 };

AtomKeySystem::AtomKeySystem(void)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::AtomKeySystem() called");
    const auto & build_in_element_map{ AtomicInfoHelper::GetElementMap() };
    const auto & build_in_remoteness_map{ AtomicInfoHelper::GetRemotenessMap() };
    const auto & build_in_branch_map{ AtomicInfoHelper::GetBranchMap() };
    for (const auto & [id0, element] : build_in_element_map)
    {
        for (const auto & [id1, remoteness] : build_in_remoteness_map)
        {
            auto id1_string{ (id1 == " ") ? "" : id1 };
            for (const auto & [id2, branch] : build_in_branch_map)
            {
                auto id2_string{ (id2 == " ") ? "" : id2 };
                auto atom_key{ GetAtomKey(element, remoteness, branch) };
                auto atom_id{ std::string{id0} + std::string{id1_string} + std::string{id2_string} };
                m_id_to_key_map.emplace(atom_id, atom_key);
                m_key_to_id_map.emplace(atom_key, atom_id);
            }
        }
    }
}

AtomKeySystem::~AtomKeySystem(void)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::~AtomKeySystem() called");
}

AtomKeySystem & AtomKeySystem::Instance(void)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::Instance() called");
    static AtomKeySystem instance;
    return instance;
}

void AtomKeySystem::RegisterAtom(const std::string & atom_id)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::RegisterAtom() called for " + atom_id);
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(atom_id) != m_id_to_key_map.end()) return;
    AtomKey new_atom_key{ m_next_dynamic_key++ };
    m_id_to_key_map[atom_id] = new_atom_key;
    m_key_to_id_map[new_atom_key] = atom_id;
}

void AtomKeySystem::RegisterAtom(const std::string & atom_id, AtomKey atom_key)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::RegisterAtom() called for " + atom_id);
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

AtomKey AtomKeySystem::GetAtomKey(Element element, Remoteness remoteness, Branch branch)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    auto key0{ static_cast<AtomKey>(element) };
    auto key1{ static_cast<AtomKey>(remoteness) };
    auto key2{ static_cast<AtomKey>(branch) };
    auto atom_key{ static_cast<AtomKey>((key2 << 16) | (key1 << 8) | key0) };
    return atom_key;
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
    return m_id_to_key_map.at(atom_id) < kDynamicBase;
}

bool AtomKeySystem::IsBuildInAtom(AtomKey atom_key) const
{
    return atom_key < kDynamicBase;
}

void AtomKeySystem::ParseAtomId(
    const std::string & atom_id,
    const std::string & element_type,
    bool is_standard_monomer,
    Element & element,
    Remoteness & remoteness,
    Branch & branch)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::ParseAtomId() called");

    element = AtomicInfoHelper::GetElementFromString(element_type);

    if (is_standard_monomer == false)
    {
        remoteness = Remoteness::UNK;
        branch = Branch::UNK;
        return;
    }

    auto atom_id_remove_element{ atom_id.substr(element_type.size()) };

    if (atom_id == "P" || atom_id_remove_element.find("P") != std::string::npos)
    {
        remoteness = Remoteness::ACID; // phosphoric acid
        std::string acid_branch_str{ atom_id_remove_element.back() };
        branch = AtomicInfoHelper::GetBranchFromString(acid_branch_str);
        return;
    }

    if (atom_id.find("'") != std::string::npos)
    {
        remoteness = Remoteness::PENTOSE; // pentose sugar (Ribose, Deoxyribose)
        auto pentose_branch_str{ atom_id_remove_element };
        pentose_branch_str.erase(
            std::remove(pentose_branch_str.begin(), pentose_branch_str.end(), '\''), pentose_branch_str.end());
        branch = AtomicInfoHelper::GetBranchFromString(pentose_branch_str);
        return;
    }

    auto remoteness_str{ StringHelper::ExtractCharAsString(atom_id_remove_element, 0) };
    if (AtomicInfoHelper::IsValidRemotenessString(remoteness_str) == true)
    {
        remoteness = AtomicInfoHelper::GetRemotenessFromString(remoteness_str);
        auto branch_str{ StringHelper::ExtractCharAsString(atom_id_remove_element, 1) };
        branch = AtomicInfoHelper::GetBranchFromString(branch_str);
        return;
    }
    else
    {
        remoteness = Remoteness::BASE; // nucleotide base
        std::string base_branch_str{ atom_id_remove_element.back() };
        branch = AtomicInfoHelper::GetBranchFromString(base_branch_str);
        return;
    }
}
