#include <gtest/gtest.h>
#include <algorithm>

#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "Logger.hpp"

class AtomicInfoHelperTest : public ::testing::Test
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

TEST_F(AtomicInfoHelperTest, BasicChecks)
{
    EXPECT_EQ(1, AtomicInfoHelper::GetAtomicNumber(Element::HYDROGEN));
    EXPECT_EQ(30, AtomicInfoHelper::GetAtomicNumber(Element::ZINC));
    EXPECT_EQ(0, AtomicInfoHelper::GetAtomicNumber(Element::UNK));
    EXPECT_EQ(0, AtomicInfoHelper::GetAtomicNumber(static_cast<Element>(999)));
    EXPECT_EQ(3U, AtomicInfoHelper::GetGroupClassCount());
    EXPECT_EQ(AtomicInfoHelper::GetStandardResidueList().size(),
              AtomicInfoHelper::GetStandardResidueCount());
    EXPECT_EQ(30U, AtomicInfoHelper::GetElementCount());
    EXPECT_EQ("element_class", AtomicInfoHelper::GetGroupClassKey(0));
    EXPECT_EQ("residue_class", AtomicInfoHelper::GetGroupClassKey(1));
    EXPECT_EQ("structure_class", AtomicInfoHelper::GetGroupClassKey(2));
}

TEST_F(AtomicInfoHelperTest, ExplicitClassKeyGetters)
{
    EXPECT_EQ("element_class", AtomicInfoHelper::GetElementClassKey());
    EXPECT_EQ("residue_class", AtomicInfoHelper::GetResidueClassKey());
    EXPECT_EQ("structure_class", AtomicInfoHelper::GetStructureClassKey());
}

TEST_F(AtomicInfoHelperTest, StandardResidueList)
{
    const auto & residue_list{ AtomicInfoHelper::GetStandardResidueList() };
    EXPECT_EQ(20u, residue_list.size());
    EXPECT_TRUE(std::find(residue_list.begin(), residue_list.end(), Residue::ALA) != residue_list.end());
}

TEST_F(AtomicInfoHelperTest, StandardElementList)
{
    const auto & element_list{ AtomicInfoHelper::GetStandardElementList() };
    EXPECT_EQ(5u, element_list.size());
    EXPECT_TRUE(std::find(element_list.begin(), element_list.end(), Element::CARBON) != element_list.end());
}

TEST_F(AtomicInfoHelperTest, StandardRemotenessList)
{
    const auto & remoteness_list{ AtomicInfoHelper::GetStandardRemotenessList() };
    EXPECT_EQ(8u, remoteness_list.size());
    EXPECT_TRUE(std::find(remoteness_list.begin(), remoteness_list.end(), Remoteness::ALPHA) != remoteness_list.end());
}

TEST_F(AtomicInfoHelperTest, StandardBranchList)
{
    const auto & branch_list{ AtomicInfoHelper::GetStandardBranchList() };
    EXPECT_EQ(4u, branch_list.size());
    EXPECT_TRUE(std::find(branch_list.begin(), branch_list.end(), Branch::ONE) != branch_list.end());
}

TEST_F(AtomicInfoHelperTest, ElementLabelMapIncludesHydrogen)
{
    const auto & label_map{ AtomicInfoHelper::GetElementLabelMap() };
    auto iter{ label_map.find(Element::HYDROGEN) };
    ASSERT_NE(iter, label_map.end());
    EXPECT_EQ("H", iter->second);
}

TEST_F(AtomicInfoHelperTest, IsStandardElement)
{
    EXPECT_TRUE(AtomicInfoHelper::IsStandardElement(Element::OXYGEN));
    EXPECT_FALSE(AtomicInfoHelper::IsStandardElement(Element::ZINC));
}

TEST_F(AtomicInfoHelperTest, IsStandardResidue)
{
    EXPECT_TRUE(AtomicInfoHelper::IsStandardResidue(Residue::ALA));
    EXPECT_FALSE(AtomicInfoHelper::IsStandardResidue(Residue::HOH));
}

TEST_F(AtomicInfoHelperTest, MapsKnownStringsToEnums)
{
    EXPECT_EQ(AtomicInfoHelper::GetResidueFromString("ALA"), Residue::ALA);
    EXPECT_EQ(AtomicInfoHelper::GetElementFromString("H"), Element::HYDROGEN);
    EXPECT_EQ(AtomicInfoHelper::GetRemotenessFromString("A"), Remoteness::ALPHA);
    EXPECT_EQ(AtomicInfoHelper::GetBranchFromString("1"), Branch::ONE);
    EXPECT_EQ(AtomicInfoHelper::GetStructureFromString("BEND"), Structure::BEND);
    EXPECT_EQ(AtomicInfoHelper::GetEntityFromString("polymer"), Entity::POLYMER);
}

TEST_F(AtomicInfoHelperTest, MapsUnknownStringsToUNK)
{
    EXPECT_EQ(AtomicInfoHelper::GetResidueFromString("XXX"), Residue::UNK);
    EXPECT_EQ(AtomicInfoHelper::GetElementFromString("Q"), Element::UNK);
    EXPECT_EQ(AtomicInfoHelper::GetRemotenessFromString("XXX"), Remoteness::UNK);
    EXPECT_EQ(AtomicInfoHelper::GetBranchFromString("XXX"), Branch::UNK);
    EXPECT_EQ(AtomicInfoHelper::GetStructureFromString("XXX"), Structure::UNK);
    EXPECT_EQ(AtomicInfoHelper::GetEntityFromString("XXX"), Entity::UNK);
}

TEST_F(AtomicInfoHelperTest, KnownResidueReturnsLabel)
{
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Residue::ALA), "ALA");
}

TEST_F(AtomicInfoHelperTest, KnownElementReturnsLabel)
{
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Element::CARBON), "C");
}

TEST_F(AtomicInfoHelperTest, ResidueUNKReturnsLabel)
{
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Residue::UNK), "UNK");
}

TEST_F(AtomicInfoHelperTest, UnknownResidueReturnsQuestionMark)
{
    EXPECT_EQ(AtomicInfoHelper::GetLabel(static_cast<Residue>(999)), "?");
}

TEST_F(AtomicInfoHelperTest, UnknownElementReturnsQuestionMark)
{
    EXPECT_EQ(AtomicInfoHelper::GetLabel(static_cast<Element>(999)), "?");
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Element::UNK), "?");
}

TEST_F(AtomicInfoHelperTest, RemotenessUNKReturnsLabel)
{
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Remoteness::UNK), "UNK");
}

TEST_F(AtomicInfoHelperTest, BranchLabel)
{
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Branch::ONE), "1");
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Branch::TERMINAL), "T");
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Branch::UNK), "UNK");
    EXPECT_EQ(AtomicInfoHelper::GetLabel(static_cast<Branch>(253)), "?");
}

TEST_F(AtomicInfoHelperTest, ElementDisplayAttributes)
{
    EXPECT_NE(1, AtomicInfoHelper::GetDisplayColor(Element::CARBON));
    EXPECT_NE(1, AtomicInfoHelper::GetDisplayMarker(Element::CARBON));
    EXPECT_EQ(1, AtomicInfoHelper::GetDisplayColor(Element::UNK));
    EXPECT_EQ(1, AtomicInfoHelper::GetDisplayMarker(Element::UNK));
    const Element invalid_element{ static_cast<Element>(999) };
    EXPECT_EQ(1, AtomicInfoHelper::GetDisplayColor(invalid_element));
    EXPECT_EQ(1, AtomicInfoHelper::GetDisplayMarker(invalid_element));
}

TEST_F(AtomicInfoHelperTest, ResidueDisplayAttributes)
{
    EXPECT_NE(1, AtomicInfoHelper::GetDisplayColor(Residue::ALA));
    EXPECT_NE(1, AtomicInfoHelper::GetDisplayMarker(Residue::ALA));
    EXPECT_EQ(1, AtomicInfoHelper::GetDisplayColor(Residue::UNK));
    EXPECT_EQ(1, AtomicInfoHelper::GetDisplayMarker(Residue::UNK));
    const Residue invalid_residue{ static_cast<Residue>(999) };
    EXPECT_EQ(1, AtomicInfoHelper::GetDisplayColor(invalid_residue));
    EXPECT_EQ(1, AtomicInfoHelper::GetDisplayMarker(invalid_residue));
}

TEST_F(AtomicInfoHelperTest, GetGroupClassKeyThrowsOnOutOfRange)
{
    EXPECT_THROW(AtomicInfoHelper::GetGroupClassKey(3), std::out_of_range);
}

TEST_F(AtomicInfoHelperTest, ConversionHelpersReturnUnkForUnknownStrings)
{
    const std::string unknown{"unknown"};
    EXPECT_EQ(Residue::UNK, AtomicInfoHelper::GetResidueFromString(unknown));
    EXPECT_EQ(Residue::UNK, AtomicInfoHelper::GetResidueFromString(unknown));
    EXPECT_EQ(Element::UNK, AtomicInfoHelper::GetElementFromString(unknown));
    EXPECT_EQ(Element::UNK, AtomicInfoHelper::GetElementFromString(unknown));
    EXPECT_EQ(Remoteness::UNK, AtomicInfoHelper::GetRemotenessFromString(unknown));
    EXPECT_EQ(Remoteness::UNK, AtomicInfoHelper::GetRemotenessFromString(unknown));
    EXPECT_EQ(Branch::UNK, AtomicInfoHelper::GetBranchFromString(unknown));
    EXPECT_EQ(Branch::UNK, AtomicInfoHelper::GetBranchFromString(unknown));
    EXPECT_EQ(Structure::UNK, AtomicInfoHelper::GetStructureFromString(unknown));
    EXPECT_EQ(Structure::UNK, AtomicInfoHelper::GetStructureFromString(unknown));
    EXPECT_EQ(Entity::UNK, AtomicInfoHelper::GetEntityFromString(unknown));
    EXPECT_EQ(Entity::UNK, AtomicInfoHelper::GetEntityFromString(unknown));
}

TEST_F(AtomicInfoHelperTest, RemotenessLabelVariants)
{
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Remoteness::NONE), "");
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Remoteness::ONE), "1");
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Remoteness::ALPHA), "#alpha");
    EXPECT_EQ(AtomicInfoHelper::GetLabel(static_cast<Remoteness>(253)), "?");
}
