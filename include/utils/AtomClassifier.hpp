#pragma once

#include <cstdint>
#include <vector>
#include <string>

enum class Element : uint16_t;
enum class Remoteness : uint8_t;
enum class Residue : uint16_t;
enum class Structure : uint8_t;
class AtomClassifier
{
    static const size_t m_main_chain_member_count{ 4 };
    static const std::vector<Element> m_main_chain_member_element_list;
    static const std::vector<Remoteness> m_main_chain_member_remoteness_list;
    static const std::vector<std::string> m_main_chain_member_label;

public:
    AtomClassifier(void);
    ~AtomClassifier();

    static bool IsMainChainMember(Element element, Remoteness remoteness, size_t & main_chain_member_id);
    static size_t GetMainChainMemberCount(void);
    static Element GetMainChainElement(size_t id);
    static Remoteness GetMainChainRemoteness(size_t id);
    static const std::string & GetMainChainElementLabel(size_t id);

    uint64_t GetMainChainElementClassGroupKey(size_t id) const;
    uint64_t GetMainChainResidueClassGroupKey(size_t id, Residue residue) const;
    uint64_t GetMainChainStructureClassGroupKey(size_t id, Structure structure, Residue residue) const;
    std::vector<uint64_t> GetMainChainResidueClassGroupKeyList(size_t id) const;
    std::vector<uint64_t> GetMainChainStructureClassGroupKeyList(size_t id, Structure structure) const;

private:


};