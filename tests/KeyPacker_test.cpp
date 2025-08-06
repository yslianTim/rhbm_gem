#include <gtest/gtest.h>
#include <unordered_set>
#include <array>

#include "KeyPacker.hpp"
#include "GlobalEnumClass.hpp"

TEST(KeyPackerElementClass, PackUnpackRoundTrip)
{
    const std::vector<Element> elements{
        Element::CARBON, Element::NITROGEN, Element::OXYGEN, Element::ZINC, Element::UNK };
    const std::vector<Remoteness> remotenesses{
        Remoteness::ALPHA, Remoteness::BETA, Remoteness::GAMMA, Remoteness::ZETA, Remoteness::UNK };

    for (auto element : elements)
    {
        for (auto remoteness : remotenesses)
        {
            for (bool flag : { false, true })
            {
                auto key{ KeyPackerElementClass::Pack(element, remoteness, flag) };
                auto [element_out, remoteness_out, flag_out]{ KeyPackerElementClass::Unpack(key) };
                EXPECT_EQ(element_out, element);
                EXPECT_EQ(remoteness_out, remoteness);
                EXPECT_EQ(flag_out, flag);
            }
        }
    }
}

TEST(KeyPackerElementClass, BoundaryValues)
{
    Element element{ Element::UNK };
    Remoteness remoteness{ Remoteness::UNK };

    auto key_false{ KeyPackerElementClass::Pack(element, remoteness, false) };
    auto expected_false{ static_cast<uint64_t>(element) | (static_cast<uint64_t>(remoteness) << 16) };
    EXPECT_EQ(key_false, expected_false);
    auto [element0, remoteness0, flag0]{ KeyPackerElementClass::Unpack(key_false) };
    EXPECT_EQ(element0, element);
    EXPECT_EQ(remoteness0, remoteness);
    EXPECT_FALSE(flag0);

    auto key_true{ KeyPackerElementClass::Pack(element, remoteness, true) };
    auto expected_true{ expected_false | (1ULL << 24) };
    EXPECT_EQ(key_true, expected_true);
    auto [element1, remoteness1, flag1]{ KeyPackerElementClass::Unpack(key_true) };
    EXPECT_EQ(element1, element);
    EXPECT_EQ(remoteness1, remoteness);
    EXPECT_TRUE(flag1);
}

TEST(KeyPackerElementClass, MinValues)
{
    Element element{ static_cast<Element>(0) };
    Remoteness remoteness{ Remoteness::NONE };

    auto key_false{ KeyPackerElementClass::Pack(element, remoteness, false) };
    EXPECT_EQ(key_false, 0u);
    auto [element0, remoteness0, flag0]{ KeyPackerElementClass::Unpack(key_false) };
    EXPECT_EQ(element0, element);
    EXPECT_EQ(remoteness0, remoteness);
    EXPECT_FALSE(flag0);

    auto key_true{ KeyPackerElementClass::Pack(element, remoteness, true) };
    EXPECT_EQ(key_true, 1ULL << 24);
    auto [element1, remoteness1, flag1]{ KeyPackerElementClass::Unpack(key_true) };
    EXPECT_EQ(element1, element);
    EXPECT_EQ(remoteness1, remoteness);
    EXPECT_TRUE(flag1);
}

TEST(KeyPackerElementClass, UniqueKeys)
{
    const std::array elements{ Element::HYDROGEN, Element::CARBON, Element::NITROGEN };
    const std::array remotenesses{ Remoteness::NONE, Remoteness::ALPHA, Remoteness::BETA };
    std::unordered_set<uint64_t> keys;
    for (auto e : elements)
    {
        for (auto r : remotenesses)
        {
            for (bool flag : { false, true })
            {
                keys.insert(KeyPackerElementClass::Pack(e, r, flag));
            }
        }
    }
    EXPECT_EQ(keys.size(), elements.size() * remotenesses.size() * 2u);
}

TEST(KeyPackerElementClass, HandlesUnknownValues)
{
    auto element{ static_cast<Element>(0x1234) };
    auto remoteness{ static_cast<Remoteness>(0xAB) };
    bool flag{ true };
    auto key{ KeyPackerElementClass::Pack(element, remoteness, flag) };
    auto [element_out, remoteness_out, flag_out]{ KeyPackerElementClass::Unpack(key) };
    EXPECT_EQ(element_out, element);
    EXPECT_EQ(remoteness_out, remoteness);
    EXPECT_EQ(flag_out, flag);
}

TEST(KeyPackerResidueClass, PackUnpackRoundTrip)
{
    const std::vector<Residue> residues{
        Residue::ALA, Residue::GLY, Residue::ZN, Residue::UNK };
    const std::vector<Element> elements{
        Element::CARBON, Element::NITROGEN, Element::ZINC, Element::UNK };
    const std::vector<Remoteness> remotenesses{
        Remoteness::ALPHA, Remoteness::BETA, Remoteness::ZETA, Remoteness::UNK };
    const std::vector<Branch> branches{
        Branch::ONE, Branch::TWO, Branch::TERMINAL, Branch::UNK };
    for (auto residue : residues)
    {
        for (auto element : elements)
        {
            for (auto remoteness : remotenesses)
            {
                for (auto branch : branches)
                {
                    for (bool flag : { false, true })
                    {
                        auto key{ KeyPackerResidueClass::Pack(
                            residue, element, remoteness, branch, flag) };
                        auto [residue_u, element_u, remoteness_u, branch_u, flag_u]{
                            KeyPackerResidueClass::Unpack(key) };
                        EXPECT_EQ(residue_u, residue);
                        EXPECT_EQ(element_u, element);
                        EXPECT_EQ(remoteness_u, remoteness);
                        EXPECT_EQ(branch_u, branch);
                        EXPECT_EQ(flag_u, flag);
                    }
                }
            }
        }
    }
}

TEST(KeyPackerResidueClass, BoundaryValues)
{
    auto residue{ Residue::UNK };
    auto element{ Element::UNK };
    auto remoteness{ Remoteness::UNK };
    auto branch{ Branch::UNK };

    auto key_false{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, false) };
    auto expected_false{ static_cast<uint64_t>(residue)
        | (static_cast<uint64_t>(element)      << 16)
        | (static_cast<uint64_t>(remoteness)   << 32)
        | (static_cast<uint64_t>(branch)       << 40) };
    EXPECT_EQ(key_false, expected_false);
    auto [residue0, element0, remoteness0, branch0, flag0]{ KeyPackerResidueClass::Unpack(key_false) };
    EXPECT_EQ(residue0, residue);
    EXPECT_EQ(element0, element);
    EXPECT_EQ(remoteness0, remoteness);
    EXPECT_EQ(branch0, branch);
    EXPECT_FALSE(flag0);

    auto key_true{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, true) };
    auto expected_true{ expected_false | (1ULL << 48) };
    EXPECT_EQ(key_true, expected_true);
    auto [residue1, element1, remoteness1, branch1, flag1]{ KeyPackerResidueClass::Unpack(key_true) };
    EXPECT_EQ(residue1, residue);
    EXPECT_EQ(element1, element);
    EXPECT_EQ(remoteness1, remoteness);
    EXPECT_EQ(branch1, branch);
    EXPECT_TRUE(flag1);
}

TEST(KeyPackerResidueClass, EdgeCases)
{
    auto residue{ Residue::UNK };
    auto element{ Element::UNK };
    auto remoteness{ Remoteness::UNK };
    auto branch{ Branch::TERMINAL };
    {
        bool flag{ true };
        auto key{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, flag) };
        EXPECT_EQ(key & 0xFFFF, static_cast<uint64_t>(residue));
        EXPECT_EQ((key >> 16) & 0xFFFF, static_cast<uint64_t>(element));
        EXPECT_EQ((key >> 32) & 0xFF, static_cast<uint64_t>(remoteness));
        EXPECT_EQ((key >> 40) & 0xFF, static_cast<uint64_t>(branch));
        EXPECT_EQ((key >> 48) & 0x1, 1u);

        auto [residue_u, element_u, remoteness_u, branch_u, flag_u]{ KeyPackerResidueClass::Unpack(key) };
        EXPECT_EQ(residue_u, residue);
        EXPECT_EQ(element_u, element);
        EXPECT_EQ(remoteness_u, remoteness);
        EXPECT_EQ(branch_u, branch);
        EXPECT_TRUE(flag_u);
    }

    {
        bool flag{ false };
        auto key{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, flag) };
        EXPECT_EQ((key >> 48) & 0x1, 0u);
        auto [residue_u, element_u, remoteness_u, branch_u, flag_u]{ KeyPackerResidueClass::Unpack(key) };
        EXPECT_FALSE(flag_u);
        EXPECT_EQ(residue_u, residue);
        EXPECT_EQ(element_u, element);
        EXPECT_EQ(remoteness_u, remoteness);
        EXPECT_EQ(branch_u, branch);
    }
}

TEST(KeyPackerResidueClass, UniqueKeys)
{
    const std::array residues{ Residue::ALA, Residue::GLY };
    const std::array elements{ Element::CARBON, Element::NITROGEN };
    const std::array remotenesses{ Remoteness::ALPHA, Remoteness::BETA };
    const std::array branches{ Branch::NONE, Branch::ONE };
    std::unordered_set<uint64_t> keys;
    for (auto res : residues)
    {
        for (auto elem : elements)
        {
            for (auto rem : remotenesses)
            {
                for (auto br : branches)
                {
                    for (bool flag : { false, true })
                    {
                        keys.insert(KeyPackerResidueClass::Pack(res, elem, rem, br, flag));
                    }
                }
            }
        }
    }

    EXPECT_EQ(keys.size(), residues.size() * elements.size() * remotenesses.size() * branches.size() * 2u);
}

TEST(KeyPackerResidueClass, HandlesUnknownValues)
{
    auto residue{ static_cast<Residue>(0x1234) };
    auto element{ static_cast<Element>(0x5678) };
    auto remoteness{ static_cast<Remoteness>(0x9A) };
    auto branch{ static_cast<Branch>(0xAB) };
    bool flag{ false };
    auto key{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, flag) };
    auto [residue_u, element_u, remoteness_u, branch_u, flag_u]{ KeyPackerResidueClass::Unpack(key) };
    EXPECT_EQ(residue_u, residue);
    EXPECT_EQ(element_u, element);
    EXPECT_EQ(remoteness_u, remoteness);
    EXPECT_EQ(branch_u, branch);
    EXPECT_EQ(flag_u, flag);
}

TEST(KeyPackerStructureClass, PackUnpackRoundTrip)
{
    const std::vector<Structure> structures{
        Structure::FREE, Structure::BEND, Structure::HELX_P, Structure::SHEET, Structure::UNK };
    const std::vector<Residue> residues{
        Residue::ALA, Residue::GLY, Residue::ZN, Residue::UNK };
    const std::vector<Element> elements{
        Element::CARBON, Element::NITROGEN, Element::ZINC, Element::UNK };
    const std::vector<Remoteness> remotenesses{
        Remoteness::ALPHA, Remoteness::BETA, Remoteness::ZETA, Remoteness::UNK };
    const std::vector<Branch> branches{
        Branch::NONE, Branch::ONE, Branch::TERMINAL, Branch::UNK };
    for (auto structure : structures)
    {
        for (auto residue : residues)
        {
            for (auto element : elements)
            {
                for (auto remoteness : remotenesses)
                {
                    for (auto branch : branches)
                    {
                        for (bool flag : { false, true })
                        {
                            auto packed{ KeyPackerStructureClass::Pack(
                                structure, residue, element, remoteness, branch, flag) };
                            auto [structure_u, residue_u, element_u, remoteness_u, branch_u, flag_u]{
                                KeyPackerStructureClass::Unpack(packed) };
                            EXPECT_EQ(structure_u, structure);
                            EXPECT_EQ(residue_u, residue);
                            EXPECT_EQ(element_u, element);
                            EXPECT_EQ(remoteness_u, remoteness);
                            EXPECT_EQ(branch_u, branch);
                            EXPECT_EQ(flag_u, flag);
                        }
                    }
                }
            }
        }
    }
}

TEST(KeyPackerStructureClass, BoundaryValues)
{
    auto structure{ Structure::UNK };
    auto residue{ Residue::UNK };
    auto element{ Element::UNK };
    auto remoteness{ Remoteness::UNK };
    auto branch{ Branch::UNK };

    auto expected_false{ static_cast<uint64_t>(structure)
        | (static_cast<uint64_t>(residue)      << 8)
        | (static_cast<uint64_t>(element)      << 24)
        | (static_cast<uint64_t>(remoteness)   << 40)
        | (static_cast<uint64_t>(branch)       << 48) };

    auto key_false{ KeyPackerStructureClass::Pack(structure, residue, element, remoteness, branch, false) };
    EXPECT_EQ(key_false, expected_false);
    auto [structure0, residue0, element0, remoteness0, branch0, flag0]{
        KeyPackerStructureClass::Unpack(key_false) };
    EXPECT_EQ(structure0, structure);
    EXPECT_EQ(residue0, residue);
    EXPECT_EQ(element0, element);
    EXPECT_EQ(remoteness0, remoteness);
    EXPECT_EQ(branch0, branch);
    EXPECT_FALSE(flag0);

    auto key_true{ KeyPackerStructureClass::Pack(structure, residue, element, remoteness, branch, true) };
    auto expected_true{ expected_false | (1ULL << 56) };
    EXPECT_EQ(key_true, expected_true);
    auto [structure1, residue1, element1, remoteness1, branch1, flag1]{
        KeyPackerStructureClass::Unpack(key_true) };
    EXPECT_EQ(structure1, structure);
    EXPECT_EQ(residue1, residue);
    EXPECT_EQ(element1, element);
    EXPECT_EQ(remoteness1, remoteness);
    EXPECT_EQ(branch1, branch);
    EXPECT_TRUE(flag1);
}

TEST(KeyPackerStructureClass, UniqueKeys)
{
    const std::array structures{ Structure::FREE, Structure::HELX_P };
    const std::array residues{ Residue::ALA, Residue::GLY };
    const std::array elements{ Element::CARBON, Element::NITROGEN };
    const std::array remotenesses{ Remoteness::ALPHA, Remoteness::BETA };
    const std::array branches{ Branch::NONE, Branch::ONE };
    std::unordered_set<uint64_t> keys;

    for (auto s : structures)
    {
        for (auto res : residues)
        {
            for (auto elem : elements)
            {
                for (auto rem : remotenesses)
                {
                    for (auto br : branches)
                    {
                        for (bool flag : { false, true })
                        {
                            keys.insert(KeyPackerStructureClass::Pack(
                                s, res, elem, rem, br, flag));
                        }
                    }
                }
            }
        }
    }

    EXPECT_EQ(keys.size(),
        structures.size() * residues.size() * elements.size() * remotenesses.size() * branches.size() * 2u);
}

TEST(KeyPackerStructureClass, EdgeCases)
{
    auto structure{Structure::UNK};
    auto residue{Residue::UNK};
    auto element{Element::UNK};
    auto remoteness{Remoteness::UNK};
    auto branch{Branch::TERMINAL};

    {
        bool flag{true};
        auto key{KeyPackerStructureClass::Pack(structure, residue, element, remoteness, branch, flag)};
        EXPECT_EQ(key & 0xFF, static_cast<uint64_t>(structure));
        EXPECT_EQ((key >> 8) & 0xFFFF, static_cast<uint64_t>(residue));
        EXPECT_EQ((key >> 24) & 0xFFFF, static_cast<uint64_t>(element));
        EXPECT_EQ((key >> 40) & 0xFF, static_cast<uint64_t>(remoteness));
        EXPECT_EQ((key >> 48) & 0xFF, static_cast<uint64_t>(branch));
        EXPECT_EQ((key >> 56) & 0x1, 1u);

        auto [structure_u, residue_u, element_u, remoteness_u, branch_u, flag_u]{ KeyPackerStructureClass::Unpack(key) };
        EXPECT_EQ(structure_u, structure);
        EXPECT_EQ(residue_u, residue);
        EXPECT_EQ(element_u, element);
        EXPECT_EQ(remoteness_u, remoteness);
        EXPECT_EQ(branch_u, branch);
        EXPECT_TRUE(flag_u);
    }

    {
        bool flag{ false };
        auto key{ KeyPackerStructureClass::Pack(structure, residue, element, remoteness, branch, flag) };
        EXPECT_EQ(key & 0xFF, static_cast<uint64_t>(structure));
        EXPECT_EQ((key >> 8) & 0xFFFF, static_cast<uint64_t>(residue));
        EXPECT_EQ((key >> 24) & 0xFFFF, static_cast<uint64_t>(element));
        EXPECT_EQ((key >> 40) & 0xFF, static_cast<uint64_t>(remoteness));
        EXPECT_EQ((key >> 48) & 0xFF, static_cast<uint64_t>(branch));
        EXPECT_EQ((key >> 56) & 0x1, 0u);

        auto [structure_u, residue_u, element_u, remoteness_u, branch_u, flag_u]{KeyPackerStructureClass::Unpack(key)};
        EXPECT_EQ(structure_u, structure);
        EXPECT_EQ(residue_u, residue);
        EXPECT_EQ(element_u, element);
        EXPECT_EQ(remoteness_u, remoteness);
        EXPECT_EQ(branch_u, branch);
        EXPECT_FALSE(flag_u);
    }
}

TEST(KeyPackerStructureClass, HandlesUnknownValues)
{
    auto structure{ static_cast<Structure>(0xCD) };
    auto residue{ static_cast<Residue>(0x1234) };
    auto element{ static_cast<Element>(0x5678) };
    auto remoteness{ static_cast<Remoteness>(0x9A) };
    auto branch{ static_cast<Branch>(0xBC) };
    bool flag{ true };

    auto key{ KeyPackerStructureClass::Pack(structure, residue, element, remoteness, branch, flag) };
    auto [structure_u, residue_u, element_u, remoteness_u, branch_u, flag_u]{
        KeyPackerStructureClass::Unpack(key) };
    EXPECT_EQ(structure_u, structure);
    EXPECT_EQ(residue_u, residue);
    EXPECT_EQ(element_u, element);
    EXPECT_EQ(remoteness_u, remoteness);
    EXPECT_EQ(branch_u, branch);
    EXPECT_EQ(flag_u, flag);
}
