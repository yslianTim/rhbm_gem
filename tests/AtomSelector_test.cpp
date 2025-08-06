#include <gtest/gtest.h>
#include <string>

#include "AtomSelector.hpp"
#include "GlobalEnumClass.hpp"
#include "Logger.hpp"

class AtomSelectorTest : public ::testing::Test
{
protected:
    LogLevel m_prev_level{};
    void SetUp(void) override
    {
        m_prev_level = Logger::GetLogLevel();
        Logger::SetLogLevel(LogLevel::Error);
    }

    void TearDown(void) override
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
    selector.PickRemotenessType("A");
    selector.VetoChainID("B");
    selector.VetoResidueType("GLY");
    selector.VetoElementType("O");
    selector.VetoRemotenessType("B");

    testing::internal::CaptureStdout();
    selector.Print();
    std::string output{ testing::internal::GetCapturedStdout() };

    EXPECT_NE(output.find("Atomic Picking List"), std::string::npos);
    EXPECT_NE(output.find("Atomic Vetoing List"), std::string::npos);
    EXPECT_NE(output.find(" - Chain set: A, "), std::string::npos);
    EXPECT_NE(output.find("ALA"), std::string::npos);
    EXPECT_NE(output.find(" - Element set: C, "), std::string::npos);
    EXPECT_NE(output.find(" - Remoteness set: #alpha, "), std::string::npos);
    EXPECT_NE(output.find("B"), std::string::npos);
    EXPECT_NE(output.find("GLY"), std::string::npos);
    EXPECT_NE(output.find("O"), std::string::npos);
    EXPECT_NE(output.find(" - Remoteness set: #beta, "), std::string::npos);
}

TEST_F(AtomSelectorTest, PrintOutputsHeadersWithEmptySets)
{
    AtomSelector selector;

    testing::internal::CaptureStdout();
    selector.Print();
    std::string output{ testing::internal::GetCapturedStdout() };

    EXPECT_NE(output.find("Atomic Picking List"), std::string::npos);
    EXPECT_NE(output.find("Atomic Vetoing List"), std::string::npos);
    EXPECT_NE(output.find(" - Chain set: \n"), std::string::npos);
    EXPECT_NE(output.find(" - Residue set: \n"), std::string::npos);
    EXPECT_NE(output.find(" - Element set: \n"), std::string::npos);
    EXPECT_NE(output.find(" - Remoteness set: \n"), std::string::npos);
}

TEST_F(AtomSelectorTest, DefaultSelectionIsTrue)
{
    AtomSelector selector;
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("B", Residue::GLY, Element::OXYGEN, Remoteness::BETA));
}

TEST_F(AtomSelectorTest, VetoChainID)
{
    AtomSelector selector;
    selector.VetoChainID("A,B");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("C", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoResidueType)
{
    AtomSelector selector;
    selector.VetoResidueType("ALA,GLY");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::GLY, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::LYS, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoElementType)
{
    AtomSelector selector;
    selector.VetoElementType("C,N");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::NITROGEN, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::OXYGEN, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoRemotenessType)
{
    AtomSelector selector;
    selector.VetoRemotenessType("A,B");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::BETA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::GAMMA));
}

TEST_F(AtomSelectorTest, PickChainID)
{
    AtomSelector selector;
    selector.PickChainID("A,B");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("C", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickResidueType)
{
    AtomSelector selector;
    selector.PickResidueType("ALA,GLY");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::GLY, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::LYS, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickChainIDAndResidueType)
{
    AtomSelector selector;
    selector.PickChainID("A");
    selector.PickResidueType("ALA");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::GLY, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickElementType)
{
    AtomSelector selector;
    selector.PickElementType("C,N");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::NITROGEN, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::OXYGEN, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickRemotenessType)
{
    AtomSelector selector;
    selector.PickRemotenessType("A,B");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::BETA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::GAMMA));
}

TEST_F(AtomSelectorTest, VetoPrecedenceOverPick)
{
    AtomSelector selector;
    selector.PickChainID("A,B");
    selector.VetoChainID("A");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoResiduePrecedenceOverPick)
{
    AtomSelector selector;
    selector.PickResidueType("ALA");
    selector.VetoResidueType("ALA");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoElementPrecedenceOverPick)
{
    AtomSelector selector;
    selector.PickElementType("C");
    selector.VetoElementType("C");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoRemotenessPrecedenceOverPick)
{
    AtomSelector selector;
    selector.PickRemotenessType("A");
    selector.VetoRemotenessType("A");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoUnknownResidueHandledAsUNK)
{
    AtomSelector selector;
    selector.VetoResidueType("XYZ");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::UNK, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickUnknownResidueHandledAsUNK)
{
    AtomSelector selector;
    selector.PickResidueType("XXX");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::UNK, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickUnknownElementHandledAsUNK)
{
    AtomSelector selector;
    selector.PickElementType("XX");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::UNK, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoUnknownElementHandledAsUNK)
{
    AtomSelector selector;
    selector.VetoElementType("XX");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::UNK, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoUnknownRemotenessHandledAsUNK)
{
    AtomSelector selector;
    selector.VetoRemotenessType("Y");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::UNK));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickUnknownRemotenessHandledAsUNK)
{
    AtomSelector selector;
    selector.PickRemotenessType("Y");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::UNK));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickChainIDHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.PickChainID("A,B,B, ,");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("C", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoChainIDHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.VetoChainID(" A,A,B,B,");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("B", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("C", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickResidueTypeHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.PickResidueType("ALA,GLY,GLY, ,");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::GLY, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::PHE, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoResidueTypeHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.VetoResidueType(" ALA,ALA,GLY,GLY,");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::GLY, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::PHE, Element::CARBON, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickElementTypeHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.PickElementType("C,N,N, ,");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::NITROGEN, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::OXYGEN, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, VetoElementTypeHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.VetoElementType(" C,C,N,N,");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::NITROGEN, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::OXYGEN, Remoteness::ALPHA));
}

TEST_F(AtomSelectorTest, PickRemotenessTypeHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.PickRemotenessType("A,B,B, ,");
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::BETA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::GAMMA));
}

TEST_F(AtomSelectorTest, VetoRemotenessTypeHandlesDuplicatesAndWhitespace)
{
    AtomSelector selector;
    selector.VetoRemotenessType(" A,A,B,B,");
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::ALPHA));
    EXPECT_FALSE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::BETA));
    EXPECT_TRUE(selector.GetSelectionFlag("A", Residue::ALA, Element::CARBON, Remoteness::GAMMA));
}
