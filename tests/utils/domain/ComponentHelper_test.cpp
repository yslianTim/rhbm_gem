#include <gtest/gtest.h>

#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

class ComponentHelperTest : public ::testing::Test
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

TEST_F(ComponentHelperTest, GetAtomCountResidue)
{
    EXPECT_EQ(5u, ComponentHelper::GetAtomCount(Residue::ALA));
    EXPECT_EQ(4u, ComponentHelper::GetAtomCount(Residue::GLY));
    EXPECT_THROW(ComponentHelper::GetAtomCount(Residue::UNK), std::out_of_range);
}

TEST_F(ComponentHelperTest, GetPartialChargeListReturnsExpected)
{
    const auto & free_list{
        ComponentHelper::GetPartialChargeList(Residue::ALA, Structure::FREE)
    };
    ASSERT_GE(free_list.size(), 1u);
    EXPECT_NEAR(0.560, free_list.front(), 1e-3);

    const auto &sheet_list{
        ComponentHelper::GetPartialChargeList(Residue::ALA, Structure::SHEET)
    };
    ASSERT_GE(sheet_list.size(), 1u);
    EXPECT_NEAR(-0.380, sheet_list.back(), 1e-3);

    const auto & helix_list{
        ComponentHelper::GetPartialChargeList(Residue::ALA, Structure::HELX_P)
    };
    ASSERT_GE(helix_list.size(), 1u);
    EXPECT_NEAR(0.559, helix_list.front(), 1e-3);

    const auto & helix_pp_list{
        ComponentHelper::GetPartialChargeList(Residue::ALA, Structure::HELX_RH_PP_P)
    };
    EXPECT_EQ(helix_list, helix_pp_list);

    EXPECT_THROW(
        ComponentHelper::GetPartialChargeList(Residue::ALA, Structure::BEND),
        std::out_of_range);
    EXPECT_THROW(
        ComponentHelper::GetPartialChargeList(Residue::ALA, Structure::TURN_P),
        std::out_of_range);
}

TEST_F(ComponentHelperTest, GetPartialChargeListAmberMatchesTable)
{
    const auto & list{ ComponentHelper::GetPartialChargeListAmber(Residue::ALA) };
    std::vector<double> expected{ 0.597, 0.034, -0.416, -0.568, -0.183 };
    ASSERT_EQ(expected.size(), list.size());
    for (std::size_t i = 0; i < expected.size(); i++)
    {
        EXPECT_NEAR(expected[i], list[i], 1e-3);
    }
}

TEST_F(ComponentHelperTest, GetPartialChargeListAmberUnknownResidueThrows)
{
    EXPECT_THROW(
        ComponentHelper::GetPartialChargeListAmber(Residue::UNK),
        std::out_of_range);
}
