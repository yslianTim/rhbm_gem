#include <gtest/gtest.h>
#include <string>

#include <rhbm_gem/utils/AtomSelector.hpp>
#include <rhbm_gem/utils/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/Logger.hpp>

class AtomSelectorTest : public ::testing::Test
{
protected:
    LogLevel m_prev_level{};
    void SetUp() override
    {
        m_prev_level = Logger::GetLogLevel();
        Logger::SetLogLevel(LogLevel::Error);
    }

    void TearDown() override
    {
        Logger::SetLogLevel(m_prev_level);
    }
};

TEST_F(AtomSelectorTest, PrintOutputsPopulatedSets)
{
    AtomSelector selector;
    selector.PickChainID("A");
    selector.PickResidueType("ALA");
    selector.PickElementType("C");
    selector.VetoChainID("B");
    selector.VetoResidueType("GLY");
    selector.VetoElementType("O");

    const std::string output{ selector.Describe() };
    Logger::SetLogLevel(LogLevel::Info);
    EXPECT_NO_THROW(selector.Print());

    EXPECT_NE(output.find("Atomic Picking List"), std::string::npos);
    EXPECT_NE(output.find("Atomic Vetoing List"), std::string::npos);
    EXPECT_NE(output.find(" - Chain set: A, "), std::string::npos);
    EXPECT_NE(output.find("ALA"), std::string::npos);
    EXPECT_NE(output.find(" - Element set: C, "), std::string::npos);
    EXPECT_NE(output.find("B"), std::string::npos);
    EXPECT_NE(output.find("GLY"), std::string::npos);
    EXPECT_NE(output.find("O"), std::string::npos);
}

TEST_F(AtomSelectorTest, PrintOutputsHeadersWithEmptySets)
{
    AtomSelector selector;

    const std::string output{ selector.Describe() };
    Logger::SetLogLevel(LogLevel::Info);
    EXPECT_NO_THROW(selector.Print());

    EXPECT_NE(output.find("Atomic Picking List"), std::string::npos);
    EXPECT_NE(output.find("Atomic Vetoing List"), std::string::npos);
    EXPECT_NE(output.find(" - Chain set: \n"), std::string::npos);
    EXPECT_NE(output.find(" - Residue set: \n"), std::string::npos);
    EXPECT_NE(output.find(" - Element set: \n"), std::string::npos);
}

TEST_F(AtomSelectorTest, DescribeReturnsPopulatedSets)
{
    AtomSelector selector;
    selector.PickChainID("A");
    selector.PickResidueType("ALA");
    selector.PickElementType("C");
    selector.VetoChainID("B");
    selector.VetoResidueType("GLY");
    selector.VetoElementType("O");

    const std::string description{ selector.Describe() };

    EXPECT_NE(description.find("Atomic Picking List"), std::string::npos);
    EXPECT_NE(description.find("Atomic Vetoing List"), std::string::npos);
    EXPECT_NE(description.find(" - Chain set: A, "), std::string::npos);
    EXPECT_NE(description.find("ALA"), std::string::npos);
    EXPECT_NE(description.find(" - Element set: C, "), std::string::npos);
    EXPECT_NE(description.find("B"), std::string::npos);
    EXPECT_NE(description.find("GLY"), std::string::npos);
    EXPECT_NE(description.find("O"), std::string::npos);
}

TEST_F(AtomSelectorTest, DescribeReturnsHeadersWithEmptySets)
{
    AtomSelector selector;

    const std::string description{ selector.Describe() };

    EXPECT_NE(description.find("Atomic Picking List"), std::string::npos);
    EXPECT_NE(description.find("Atomic Vetoing List"), std::string::npos);
    EXPECT_NE(description.find(" - Chain set: \n"), std::string::npos);
    EXPECT_NE(description.find(" - Residue set: \n"), std::string::npos);
    EXPECT_NE(description.find(" - Element set: \n"), std::string::npos);
}

TEST_F(AtomSelectorTest, DefaultSelectionIsTrue)
{
    AtomSelector selector;
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("B", Residue::GLY, Element::OXYGEN));
}

TEST_F(AtomSelectorTest, VetoChainID)
{
    AtomSelector selector;
    selector.VetoChainID("A,B");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("C", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, VetoResidueType)
{
    AtomSelector selector;
    selector.VetoResidueType("ALA,GLY");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::GLY, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::LYS, Element::CARBON));
}

TEST_F(AtomSelectorTest, VetoElementType)
{
    AtomSelector selector;
    selector.VetoElementType("C,N");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::NITROGEN));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::OXYGEN));
}

TEST_F(AtomSelectorTest, PickChainID)
{
    AtomSelector selector;
    selector.PickChainID("A,B");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("C", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, PickResidueType)
{
    AtomSelector selector;
    selector.PickResidueType("ALA,GLY");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::GLY, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::LYS, Element::CARBON));
}

TEST_F(AtomSelectorTest, PickChainIDAndResidueType)
{
    AtomSelector selector;
    selector.PickChainID("A");
    selector.PickResidueType("ALA");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::GLY, Element::CARBON));
}

TEST_F(AtomSelectorTest, PickElementType)
{
    AtomSelector selector;
    selector.PickElementType("C,N");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::NITROGEN));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::OXYGEN));
}

TEST_F(AtomSelectorTest, VetoPrecedenceOverPick)
{
    AtomSelector selector;
    selector.PickChainID("A,B");
    selector.VetoChainID("A");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, VetoResiduePrecedenceOverPick)
{
    AtomSelector selector;
    selector.PickResidueType("ALA");
    selector.VetoResidueType("ALA");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, VetoElementPrecedenceOverPick)
{
    AtomSelector selector;
    selector.PickElementType("C");
    selector.VetoElementType("C");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, VetoUnknownResidueHandledAsUNK)
{
    AtomSelector selector;
    selector.VetoResidueType("XYZ");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::UNK, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, PickUnknownResidueHandledAsUNK)
{
    AtomSelector selector;
    selector.PickResidueType("XXX");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::UNK, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, PickUnknownElementHandledAsUNK)
{
    AtomSelector selector;
    selector.PickElementType("XX");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::UNK));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, VetoUnknownElementHandledAsUNK)
{
    AtomSelector selector;
    selector.VetoElementType("XX");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::UNK));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, PickChainIDHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.PickChainID("A,B,B, ,");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("C", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, VetoChainIDHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.VetoChainID(" A,A,B,B,");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("C", Residue::ALA, Element::CARBON));
}

TEST_F(AtomSelectorTest, PickResidueTypeHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.PickResidueType("ALA,GLY,GLY, ,");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::GLY, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::PHE, Element::CARBON));
}

TEST_F(AtomSelectorTest, VetoResidueTypeHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.VetoResidueType(" ALA,ALA,GLY,GLY,");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::GLY, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::PHE, Element::CARBON));
}

TEST_F(AtomSelectorTest, PickElementTypeHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.PickElementType("C,N,N, ,");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::NITROGEN));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::OXYGEN));
}

TEST_F(AtomSelectorTest, VetoElementTypeHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.VetoElementType(" C,C,N,N,");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::NITROGEN));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::OXYGEN));
}
