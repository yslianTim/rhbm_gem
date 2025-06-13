#pragma once

#include <cstdint>
#include <vector>
#include <string>

enum class Element : uint16_t;
enum class Remoteness : uint8_t;
enum class Residue : uint16_t;
enum class Structure : uint8_t;

class AtomObject;
class AtomClassifier
{
    static const size_t m_main_chain_member_count{ 4 };
    static const std::vector<short> m_main_chain_member_color_list;
    static const std::vector<short> m_main_chain_member_solid_marker_list;
    static const std::vector<short> m_main_chain_member_open_marker_list;
    static const std::vector<Element> m_main_chain_member_element_list;
    static const std::vector<Remoteness> m_main_chain_member_remoteness_list;
    static const std::vector<std::string> m_main_chain_member_label_list;
    static const std::vector<std::string> m_main_chain_member_title_list;

public:
    AtomClassifier(void);
    ~AtomClassifier();

    static bool IsMainChainMember(Element element, Remoteness remoteness, size_t & main_chain_member_id);
    static size_t GetMainChainMemberCount(void);
    static Element GetMainChainElement(size_t id);
    static Remoteness GetMainChainRemoteness(size_t id);
    static short GetMainChainElementColor(size_t id);
    static short GetMainChainElementSolidMarker(size_t id);
    static short GetMainChainElementOpenMarker(size_t id);
    static const std::string & GetMainChainElementLabel(size_t id);
    static const std::string & GetMainChainElementTitle(size_t id);

    static uint64_t GetGroupKeyInClass(const AtomObject * atom_object, const std::string & class_key);

    uint64_t GetMainChainElementClassGroupKey(size_t id) const;
    uint64_t GetMainChainResidueClassGroupKey(size_t id, Residue residue) const;
    uint64_t GetMainChainStructureClassGroupKey(size_t id, Structure structure, Residue residue) const;
    std::vector<uint64_t> GetMainChainResidueClassGroupKeyList(size_t id) const;
    std::vector<uint64_t> GetMainChainStructureClassGroupKeyList(size_t id, Structure structure) const;

private:
    static bool IsValidMainChainMemberID(size_t id);

};