#include <gtest/gtest.h>

#include <set>

#include "OptionEnumTraits.hpp"

namespace rg = rhbm_gem;

namespace {

template <typename EnumType>
std::set<int> BuildValueSetFromCLI()
{
    std::set<int> values;
    for (const auto & [token, value] : rg::BuildEnumCLIMap<EnumType>())
    {
        static_cast<void>(token);
        values.insert(static_cast<int>(value));
    }
    return values;
}

template <typename EnumType>
std::set<int> BuildValueSetFromBindings()
{
    std::set<int> values;
    for (const auto & entry : rg::EnumOptionTraits<EnumType>::kBindingEntries)
    {
        values.insert(static_cast<int>(entry.value));
    }
    return values;
}

} // namespace

TEST(EnumOptionTraitsTest, PainterTypeTraitsStayInSync)
{
    const auto cli_map{ rg::BuildEnumCLIMap<rg::PainterType>() };
    EXPECT_EQ(cli_map.at("gaus"), rg::PainterType::GAUS);
    EXPECT_EQ(cli_map.at("atom"), rg::PainterType::ATOM);
    EXPECT_EQ(BuildValueSetFromCLI<rg::PainterType>(), BuildValueSetFromBindings<rg::PainterType>());
}

TEST(EnumOptionTraitsTest, PrinterTypeTraitsStayInSync)
{
    const auto cli_map{ rg::BuildEnumCLIMap<rg::PrinterType>() };
    EXPECT_EQ(cli_map.at("atom_out"), rg::PrinterType::ATOM_OUTLIER);
    EXPECT_EQ(
        BuildValueSetFromCLI<rg::PrinterType>(),
        BuildValueSetFromBindings<rg::PrinterType>());
}

TEST(EnumOptionTraitsTest, PotentialModelTraitsStayInSync)
{
    const auto cli_map{ rg::BuildEnumCLIMap<rg::PotentialModel>() };
    EXPECT_EQ(cli_map.at("five"), rg::PotentialModel::FIVE_GAUS_CHARGE);
    EXPECT_EQ(
        BuildValueSetFromCLI<rg::PotentialModel>(),
        BuildValueSetFromBindings<rg::PotentialModel>());
}

TEST(EnumOptionTraitsTest, PartialChargeTraitsStayInSync)
{
    const auto cli_map{ rg::BuildEnumCLIMap<rg::PartialCharge>() };
    EXPECT_EQ(cli_map.at("amber"), rg::PartialCharge::AMBER);
    EXPECT_EQ(
        BuildValueSetFromCLI<rg::PartialCharge>(),
        BuildValueSetFromBindings<rg::PartialCharge>());
}
