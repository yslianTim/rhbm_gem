#include <gtest/gtest.h>
#include <algorithm>
#include <utility>

#include <rhbm_gem/utils/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/Logger.hpp>

class ChemicalDataHelperTest : public ::testing::Test
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

TEST_F(ChemicalDataHelperTest, BasicChecks)
{
    EXPECT_EQ(3U, ChemicalDataHelper::GetGroupAtomClassCount());
    EXPECT_EQ(ChemicalDataHelper::GetStandardAminoAcidList().size(),
              ChemicalDataHelper::GetStandardAminoAcidCount());
    EXPECT_EQ(90U, ChemicalDataHelper::GetElementCount());
}

using ElementAtomicPair = std::pair<Element, int>;

class ChemicalDataHelperAtomicNumberTest : public ChemicalDataHelperTest,
                                         public ::testing::WithParamInterface<ElementAtomicPair>
{
};

TEST_P(ChemicalDataHelperAtomicNumberTest, ReturnsExpectedAtomicNumber)
{
    const auto [element, atomic_number]{ GetParam() };
    EXPECT_EQ(atomic_number, ChemicalDataHelper::GetAtomicNumber(element));
}

INSTANTIATE_TEST_SUITE_P(AtomicNumberPairs, ChemicalDataHelperAtomicNumberTest,
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

TEST_F(ChemicalDataHelperTest, InvalidElementReturnsZero)
{
    EXPECT_EQ(0, ChemicalDataHelper::GetAtomicNumber(static_cast<Element>(0)));
    EXPECT_EQ(0, ChemicalDataHelper::GetAtomicNumber(Element::UNK));
}

TEST_F(ChemicalDataHelperTest, ExplicitClassKeyGetters)
{
    EXPECT_EQ("simple_atom_class", ChemicalDataHelper::GetSimpleAtomClassKey());
    EXPECT_EQ("component_atom_class", ChemicalDataHelper::GetComponentAtomClassKey());
    EXPECT_EQ("structure_atom_class", ChemicalDataHelper::GetStructureAtomClassKey());
}

TEST_F(ChemicalDataHelperTest, StandardResidueList)
{
    const auto & residue_list{ ChemicalDataHelper::GetStandardAminoAcidList() };
    EXPECT_EQ(20u, residue_list.size());
    EXPECT_TRUE(std::find(residue_list.begin(), residue_list.end(), Residue::ALA) != residue_list.end());
}

TEST_F(ChemicalDataHelperTest, StandardElementList)
{
    const auto & element_list{ ChemicalDataHelper::GetStandardElementList() };
    EXPECT_EQ(5u, element_list.size());
    EXPECT_TRUE(std::find(element_list.begin(), element_list.end(), Element::CARBON) != element_list.end());
}

TEST_F(ChemicalDataHelperTest, ElementLabelMapIncludesHydrogen)
{
    const auto & label_map{ ChemicalDataHelper::GetElementLabelMap() };
    auto iter{ label_map.find(Element::HYDROGEN) };
    ASSERT_NE(iter, label_map.end());
    EXPECT_EQ("H", iter->second);
}

TEST_F(ChemicalDataHelperTest, IsStandardElement)
{
    for (auto element : ChemicalDataHelper::GetStandardElementList())
    {
        EXPECT_TRUE(ChemicalDataHelper::IsStandardElement(element)) << static_cast<int>(element);
    }
    EXPECT_FALSE(ChemicalDataHelper::IsStandardElement(Element::ZINC));
    EXPECT_FALSE(ChemicalDataHelper::IsStandardElement(Element::PHOSPHORUS));
}

TEST_F(ChemicalDataHelperTest, IsStandardResidue)
{
    for (auto residue : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        EXPECT_TRUE(ChemicalDataHelper::IsStandardAminoAcid(residue))
            << static_cast<int>(residue);
    }
    EXPECT_FALSE(ChemicalDataHelper::IsStandardAminoAcid(Residue::UNK));
}

TEST_F(ChemicalDataHelperTest, MapsKnownStringsToEnums)
{
    EXPECT_EQ(ChemicalDataHelper::GetResidueFromString("ALA"), Residue::ALA);
    EXPECT_EQ(ChemicalDataHelper::GetElementFromString("H"), Element::HYDROGEN);
    EXPECT_EQ(ChemicalDataHelper::GetStructureFromString("BEND"), Structure::BEND);
    EXPECT_EQ(ChemicalDataHelper::GetEntityFromString("polymer"), Entity::POLYMER);
}

TEST_F(ChemicalDataHelperTest, MapsUnknownStringsToUNK)
{
    EXPECT_EQ(ChemicalDataHelper::GetResidueFromString("XXX"), Residue::UNK);
    EXPECT_EQ(ChemicalDataHelper::GetElementFromString("Q"), Element::UNK);
    EXPECT_EQ(ChemicalDataHelper::GetStructureFromString("XXX"), Structure::UNK);
    EXPECT_EQ(ChemicalDataHelper::GetEntityFromString("XXX"), Entity::UNK);
}

TEST_F(ChemicalDataHelperTest, KnownResidueReturnsLabel)
{
    EXPECT_EQ(ChemicalDataHelper::GetLabel(Residue::ALA), "ALA");
}

TEST_F(ChemicalDataHelperTest, KnownElementReturnsLabel)
{
    EXPECT_EQ(ChemicalDataHelper::GetLabel(Element::CARBON), "C");
}

TEST_F(ChemicalDataHelperTest, ResidueUNKReturnsLabel)
{
    EXPECT_EQ(ChemicalDataHelper::GetLabel(Residue::UNK), "UNK");
}

TEST_F(ChemicalDataHelperTest, UnknownResidueReturnsQuestionMark)
{
    EXPECT_EQ(ChemicalDataHelper::GetLabel(static_cast<Residue>(255)), "?");
}

TEST_F(ChemicalDataHelperTest, UnknownElementReturnsQuestionMark)
{
    EXPECT_EQ(ChemicalDataHelper::GetLabel(static_cast<Element>(255)), "?");
    EXPECT_EQ(ChemicalDataHelper::GetLabel(Element::UNK), "UNK");
}

TEST_F(ChemicalDataHelperTest, ElementDisplayAttributes)
{
    EXPECT_NE(1, ChemicalDataHelper::GetDisplayColor(Element::CARBON));
    EXPECT_NE(1, ChemicalDataHelper::GetDisplayMarker(Element::CARBON));
    EXPECT_EQ(1, ChemicalDataHelper::GetDisplayColor(Element::UNK));
    EXPECT_EQ(1, ChemicalDataHelper::GetDisplayMarker(Element::UNK));
}

TEST_F(ChemicalDataHelperTest, ResidueDisplayAttributes)
{
    EXPECT_NE(1, ChemicalDataHelper::GetDisplayColor(Residue::ALA));
    EXPECT_NE(1, ChemicalDataHelper::GetDisplayMarker(Residue::ALA));
    EXPECT_EQ(1, ChemicalDataHelper::GetDisplayColor(Residue::UNK));
    EXPECT_EQ(1, ChemicalDataHelper::GetDisplayMarker(Residue::UNK));
}

TEST_F(ChemicalDataHelperTest, GetGroupAtomClassKeyThrowsOnOutOfRange)
{
    EXPECT_THROW(ChemicalDataHelper::GetGroupAtomClassKey(3), std::out_of_range);
}

TEST_F(ChemicalDataHelperTest, ConversionHelpersReturnUnkForUnknownStrings)
{
    const std::string unknown{"unknown"};
    EXPECT_EQ(Residue::UNK, ChemicalDataHelper::GetResidueFromString(unknown));
    EXPECT_EQ(Residue::UNK, ChemicalDataHelper::GetResidueFromString(unknown));
    EXPECT_EQ(Element::UNK, ChemicalDataHelper::GetElementFromString(unknown));
    EXPECT_EQ(Element::UNK, ChemicalDataHelper::GetElementFromString(unknown));
    EXPECT_EQ(Structure::UNK, ChemicalDataHelper::GetStructureFromString(unknown));
    EXPECT_EQ(Structure::UNK, ChemicalDataHelper::GetStructureFromString(unknown));
    EXPECT_EQ(Entity::UNK, ChemicalDataHelper::GetEntityFromString(unknown));
    EXPECT_EQ(Entity::UNK, ChemicalDataHelper::GetEntityFromString(unknown));
}
