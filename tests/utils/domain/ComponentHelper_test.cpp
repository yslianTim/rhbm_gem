#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <future>

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

TEST_F(ComponentHelperTest, LookupIsIndependentOfStructureAndTableCallOrder)
{
    std::array<Structure, 4> structures{
        Structure::FREE, Structure::HELX_P, Structure::SHEET, Structure::BEND };
    for (int pass = 0; pass < 2; ++pass)
    {
        for (const auto structure : structures)
        {
            EXPECT_DOUBLE_EQ(ComponentHelper::GetPartialCharge(
                Residue::ALA, Spot::O, structure, true), -0.568);
            const auto result{ ComponentHelper::LookupPartialCharge(
                Residue::ALA, Spot::O, structure) };
            if (structure == Structure::BEND)
            {
                EXPECT_EQ(result.status, ChargeLookupStatus::UnsupportedStructure);
                EXPECT_FALSE(result.charge);
                EXPECT_FALSE(result.table);
                EXPECT_DOUBLE_EQ(ComponentHelper::GetPartialCharge(
                    Residue::ALA, Spot::O, structure), 0.0);
            }
            else
            {
                ASSERT_TRUE(result.charge);
                EXPECT_EQ(result.status, ChargeLookupStatus::Found);
                EXPECT_DOUBLE_EQ(*result.charge,
                    ComponentHelper::GetPartialChargeList(Residue::ALA, structure).at(3));
            }
        }
        std::reverse(structures.begin(), structures.end());
    }
}

TEST_F(ComponentHelperTest, EveryChargeTableMatchesItsSpotList)
{
    for (int id = static_cast<int>(Residue::ALA); id <= static_cast<int>(Residue::VAL); ++id)
    {
        const auto residue{ static_cast<Residue>(id) };
        const auto & spots{ ComponentHelper::GetSpotList(residue) };
        for (const auto structure : { Structure::FREE, Structure::HELX_RH_PP_P, Structure::SHEET })
        {
            for (const bool amber : { false, true })
            {
                const auto & charges{ amber ? ComponentHelper::GetPartialChargeListAmber(residue)
                    : ComponentHelper::GetPartialChargeList(residue, structure) };
                ASSERT_EQ(spots.size(), charges.size());
                for (size_t i = 0; i < spots.size(); ++i)
                {
                    const auto result{ ComponentHelper::LookupPartialCharge(residue, spots[i], structure, amber) };
                    ASSERT_EQ(result.status, ChargeLookupStatus::Found);
                    ASSERT_TRUE(result.charge);
                    EXPECT_DOUBLE_EQ(*result.charge, charges[i]);
                    EXPECT_EQ(result.table, amber ? ChargeTable::Amber95
                        : structure == Structure::FREE ? ChargeTable::Buried
                        : structure == Structure::SHEET ? ChargeTable::Sheet : ChargeTable::Helix);
                }
            }
        }
    }
}

TEST_F(ComponentHelperTest, LookupDistinguishesZeroFromMissingCharge)
{
    const auto zero{ ComponentHelper::LookupPartialCharge(Residue::VAL, Spot::CB, Structure::SHEET) };
    ASSERT_TRUE(zero.charge);
    EXPECT_DOUBLE_EQ(*zero.charge, 0.0);
    EXPECT_EQ(zero.status, ChargeLookupStatus::Found);
    const auto unknown{ ComponentHelper::LookupPartialCharge(Residue::UNK, Spot::O, Structure::FREE) };
    EXPECT_FALSE(unknown.charge);
    EXPECT_EQ(unknown.table, ChargeTable::Buried);
    EXPECT_EQ(unknown.status, ChargeLookupStatus::UnsupportedResidue);
    const auto missing{ ComponentHelper::LookupPartialCharge(Residue::ALA, Spot::UNK, Structure::FREE) };
    EXPECT_FALSE(missing.charge);
    EXPECT_EQ(missing.status, ChargeLookupStatus::UnsupportedSpot);
}

TEST_F(ComponentHelperTest, ConcurrentLookupsDoNotShareMutableState)
{
    std::vector<std::future<bool>> workers;
    for (const auto structure : { Structure::FREE, Structure::HELX_P, Structure::SHEET })
    {
        workers.emplace_back(std::async(std::launch::async, [structure]
        {
            for (int repeat = 0; repeat < 1000; ++repeat)
            {
                const bool amber{ repeat % 2 == 0 };
                const auto expected{ amber ? -0.568
                    : ComponentHelper::GetPartialChargeList(Residue::ALA, structure).at(3) };
                if (ComponentHelper::GetPartialCharge(Residue::ALA, Spot::O, structure, amber) != expected)
                    return false;
            }
            return true;
        }));
    }
    for (auto & worker : workers) EXPECT_TRUE(worker.get());
}
