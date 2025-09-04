#include <gtest/gtest.h>
#include <algorithm>
#include <utility>

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
    EXPECT_EQ(3U, AtomicInfoHelper::GetGroupClassCount());
    EXPECT_EQ(AtomicInfoHelper::GetStandardResidueList().size(),
              AtomicInfoHelper::GetStandardResidueCount());
    EXPECT_EQ(90U, AtomicInfoHelper::GetElementCount());
    EXPECT_EQ("element_class", AtomicInfoHelper::GetGroupClassKey(0));
    EXPECT_EQ("residue_class", AtomicInfoHelper::GetGroupClassKey(1));
    EXPECT_EQ("structure_class", AtomicInfoHelper::GetGroupClassKey(2));
}

using ElementAtomicPair = std::pair<Element, int>;

class AtomicInfoHelperAtomicNumberTest : public AtomicInfoHelperTest,
                                         public ::testing::WithParamInterface<ElementAtomicPair>
{
};

TEST_P(AtomicInfoHelperAtomicNumberTest, ReturnsExpectedAtomicNumber)
{
    const auto [element, atomic_number]{ GetParam() };
    EXPECT_EQ(atomic_number, AtomicInfoHelper::GetAtomicNumber(element));
}

INSTANTIATE_TEST_SUITE_P(AtomicNumberPairs, AtomicInfoHelperAtomicNumberTest,
    ::testing::Values(
        ElementAtomicPair{ Element::HYDROGEN,   1 },
        ElementAtomicPair{ Element::HELIUM,     2 },
        ElementAtomicPair{ Element::LITHIUM,    3 },
        ElementAtomicPair{ Element::BERYLLIUM,  4 },
        ElementAtomicPair{ Element::BORON,      5 },
        ElementAtomicPair{ Element::CARBON,     6 },
        ElementAtomicPair{ Element::NITROGEN,   7 },
        ElementAtomicPair{ Element::OXYGEN,     8 },
        ElementAtomicPair{ Element::FLUORINE,   9 },
        ElementAtomicPair{ Element::NEON,      10 },
        ElementAtomicPair{ Element::SODIUM,    11 },
        ElementAtomicPair{ Element::MAGNESIUM, 12 },
        ElementAtomicPair{ Element::ALUMINUM,  13 },
        ElementAtomicPair{ Element::SILICON,   14 },
        ElementAtomicPair{ Element::PHOSPHORUS,15 },
        ElementAtomicPair{ Element::SULFUR,    16 },
        ElementAtomicPair{ Element::CHLORINE,  17 },
        ElementAtomicPair{ Element::ARGON,     18 },
        ElementAtomicPair{ Element::POTASSIUM, 19 },
        ElementAtomicPair{ Element::CALCIUM,   20 },
        ElementAtomicPair{ Element::SCANDIUM,  21 },
        ElementAtomicPair{ Element::TITANIUM,  22 },
        ElementAtomicPair{ Element::VANADIUM,  23 },
        ElementAtomicPair{ Element::CHROMIUM,  24 },
        ElementAtomicPair{ Element::MANGANESE, 25 },
        ElementAtomicPair{ Element::IRON,      26 },
        ElementAtomicPair{ Element::COBALT,    27 },
        ElementAtomicPair{ Element::NICKEL,    28 },
        ElementAtomicPair{ Element::COPPER,    29 },
        ElementAtomicPair{ Element::ZINC,      30 }
    )
);

TEST_F(AtomicInfoHelperTest, InvalidElementReturnsZero)
{
    EXPECT_EQ(0, AtomicInfoHelper::GetAtomicNumber(static_cast<Element>(255)));
    EXPECT_EQ(0, AtomicInfoHelper::GetAtomicNumber(Element::UNK));
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
    for (auto element : AtomicInfoHelper::GetStandardElementList())
    {
        EXPECT_TRUE(AtomicInfoHelper::IsStandardElement(element)) << static_cast<int>(element);
    }
    EXPECT_FALSE(AtomicInfoHelper::IsStandardElement(Element::ZINC));
    EXPECT_FALSE(AtomicInfoHelper::IsStandardElement(Element::PHOSPHORUS));
}

TEST_F(AtomicInfoHelperTest, IsStandardResidue)
{
    for (auto residue : AtomicInfoHelper::GetStandardResidueList())
    {
        EXPECT_TRUE(AtomicInfoHelper::IsStandardResidue(residue))
            << static_cast<int>(residue);
    }
    EXPECT_FALSE(AtomicInfoHelper::IsStandardResidue(Residue::UNK));
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
    EXPECT_EQ(AtomicInfoHelper::GetLabel(static_cast<Residue>(255)), "?");
}

TEST_F(AtomicInfoHelperTest, UnknownElementReturnsQuestionMark)
{
    EXPECT_EQ(AtomicInfoHelper::GetLabel(static_cast<Element>(255)), "?");
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Element::UNK), "UNK");
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
}

TEST_F(AtomicInfoHelperTest, ResidueDisplayAttributes)
{
    EXPECT_NE(1, AtomicInfoHelper::GetDisplayColor(Residue::ALA));
    EXPECT_NE(1, AtomicInfoHelper::GetDisplayMarker(Residue::ALA));
    EXPECT_EQ(1, AtomicInfoHelper::GetDisplayColor(Residue::UNK));
    EXPECT_EQ(1, AtomicInfoHelper::GetDisplayMarker(Residue::UNK));
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
    EXPECT_EQ(AtomicInfoHelper::GetLabel(Remoteness::ALPHA), "#alpha");
    EXPECT_EQ(AtomicInfoHelper::GetLabel(static_cast<Remoteness>(253)), "?");
}
