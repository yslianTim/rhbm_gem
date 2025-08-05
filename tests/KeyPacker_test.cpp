#include <gtest/gtest.h>

#include "KeyPacker.hpp"
#include "GlobalEnumClass.hpp"

TEST(KeyPackerElementClass, PackUnpackRoundTrip)
{
    Element element{ Element::CARBON };
    Remoteness remoteness{ Remoteness::ALPHA };
    for (bool flag : { false, true })
    {
        auto key{ KeyPackerElementClass::Pack(element, remoteness, flag) };
        auto [element_out, remoteness_out, flag_out]{ KeyPackerElementClass::Unpack(key) };
        EXPECT_EQ(element_out, element);
        EXPECT_EQ(remoteness_out, remoteness);
        if (flag)
        {
            EXPECT_TRUE(flag_out);
        }
        else
        {
            EXPECT_FALSE(flag_out);
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

TEST(KeyPackerResidueClass, PackUnpackRoundTrip)
{
    auto residue{ Residue::ALA };
    auto element{ Element::CARBON };
    auto remoteness{ Remoteness::ALPHA };
    auto branch{ Branch::ONE };
    bool flag{ false };
    auto key{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, flag) };
    auto [residue_u, element_u, remoteness_u, branch_u, flag_u]{ KeyPackerResidueClass::Unpack(key) };
    EXPECT_EQ(residue_u, residue);
    EXPECT_EQ(element_u, element);
    EXPECT_EQ(remoteness_u, remoteness);
    EXPECT_EQ(branch_u, branch);
    EXPECT_EQ(flag_u, flag);
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

TEST(KeyPackerStructureClass, PackUnpackRoundTrip)
{
    struct TestCase
    {
        Structure structure;
        Residue residue;
        Element element;
        Remoteness remoteness;
        Branch branch;
        bool flag;
    };

    const std::vector<TestCase> cases{
        {Structure::FREE,    Residue::ALA, Element::CARBON,   Remoteness::ALPHA, Branch::NONE,   false},
        {Structure::HELX_P,  Residue::GLY, Element::NITROGEN, Remoteness::BETA,  Branch::ONE,    true},
        {Structure::SHEET,   Residue::SER, Element::OXYGEN,   Remoteness::GAMMA, Branch::TWO,    false},
        {Structure::STRN,    Residue::VAL, Element::SULFUR,   Remoteness::DELTA, Branch::THREE,  true},
        // Boundary case with highest enum values to ensure proper bit allocation
        {Structure::TURN_P, Residue::UNK, Element::UNK, Remoteness::UNK, Branch::UNK, true}
    };

    for (const auto & c : cases)
    {
        auto packed{ KeyPackerStructureClass::Pack(
            c.structure, c.residue, c.element, c.remoteness, c.branch, c.flag) };
        auto [structure, residue, element, remoteness, branch, flag]{ KeyPackerStructureClass::Unpack(packed) };
        EXPECT_EQ(c.structure,   structure);
        EXPECT_EQ(c.residue,     residue);
        EXPECT_EQ(c.element,     element);
        EXPECT_EQ(c.remoteness,  remoteness);
        EXPECT_EQ(c.branch,      branch);
        EXPECT_EQ(c.flag,        flag);
    }
}
