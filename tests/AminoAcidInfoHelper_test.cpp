#include <gtest/gtest.h>

#include "AminoAcidInfoHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "Logger.hpp"

class AminoAcidInfoHelperTest : public ::testing::Test
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

TEST_F(AminoAcidInfoHelperTest, GetAtomCountResidue)
{
    EXPECT_EQ(5u, AminoAcidInfoHelper::GetAtomCount(Residue::ALA));
    EXPECT_EQ(4u, AminoAcidInfoHelper::GetAtomCount(Residue::GLY));
    EXPECT_THROW(AminoAcidInfoHelper::GetAtomCount(Residue::UNK), std::out_of_range);
}

TEST_F(AminoAcidInfoHelperTest, GetAtomCountInt)
{
    EXPECT_EQ(5u, AminoAcidInfoHelper::GetAtomCount(static_cast<int>(Residue::ALA)));
    EXPECT_EQ(4u, AminoAcidInfoHelper::GetAtomCount(static_cast<int>(Residue::GLY)));
    EXPECT_THROW(
        AminoAcidInfoHelper::GetAtomCount(static_cast<int>(Residue::UNK)),
        std::out_of_range);
}

TEST_F(AminoAcidInfoHelperTest, GetPartialChargeListReturnsExpected)
{
    const auto & free_list{ AminoAcidInfoHelper::GetPartialChargeList(Residue::ALA, Structure::FREE) };
    ASSERT_GE(free_list.size(), 1u);
    EXPECT_NEAR(0.560, free_list.front(), 1e-3);

    const auto &sheet_list{ AminoAcidInfoHelper::GetPartialChargeList(Residue::ALA, Structure::SHEET) };
    ASSERT_GE(sheet_list.size(), 1u);
    EXPECT_NEAR(-0.380, sheet_list.back(), 1e-3);

    const auto & helix_list{ AminoAcidInfoHelper::GetPartialChargeList(Residue::ALA, Structure::HELX_P) };
    ASSERT_GE(helix_list.size(), 1u);
    EXPECT_NEAR(0.559, helix_list.front(), 1e-3);

    EXPECT_THROW(
        AminoAcidInfoHelper::GetPartialChargeList(Residue::ALA, Structure::BEND),
        std::out_of_range);
    EXPECT_THROW(
        AminoAcidInfoHelper::GetPartialChargeList(Residue::ALA, Structure::TURN_P),
        std::out_of_range);
}

TEST_F(AminoAcidInfoHelperTest, GetPartialChargeListAmberMatchesTable)
{
    const auto & list{ AminoAcidInfoHelper::GetPartialChargeListAmber(Residue::ALA) };
    std::vector<double> expected{ 0.597, 0.034, -0.416, -0.568, -0.183 };
    ASSERT_EQ(expected.size(), list.size());
    for (std::size_t i = 0; i < expected.size(); i++)
    {
        EXPECT_NEAR(expected[i], list[i], 1e-3);
    }
}

TEST_F(AminoAcidInfoHelperTest, GetPartialChargeValidAndCached)
{
    const double q1{
        AminoAcidInfoHelper::GetPartialCharge(
            Residue::ALA, Element::CARBON, Remoteness::BETA, Branch::NONE, Structure::FREE)
    };
    EXPECT_NEAR(-0.392, q1, 1e-3);
    const double q2{
        AminoAcidInfoHelper::GetPartialCharge(
        Residue::ALA, Element::CARBON, Remoteness::BETA, Branch::NONE, Structure::FREE)};
    EXPECT_NEAR(-0.392, q2, 1e-3);
}
