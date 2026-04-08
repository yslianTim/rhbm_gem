#include <gtest/gtest.h>

#include <array>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "command/detail/CommandEnumMetadata.hpp"
#include "support/PublicHeaderSurfaceTestSupport.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/core/command/CommandEnums.hpp>

namespace rg = rhbm_gem;

namespace {

template <typename EnumType>
std::set<int> BuildValueSetFromCLI()
{
    std::set<int> values;
    for (const auto & [token, value] : rg::internal::BuildCommandEnumCliMap<EnumType>())
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
    const auto binding_entries{ rg::internal::GetCommandEnumBindingEntries<EnumType>() };
    for (const auto & entry : binding_entries)
    {
        values.insert(static_cast<int>(entry.value));
    }
    return values;
}

template <typename EnumType>
struct EnumMappingExpectation
{
    std::string_view token;
    EnumType value;
};

template <typename EnumType>
struct EnumMappingTraits;

template <>
struct EnumMappingTraits<rg::PainterType>
{
    static constexpr std::string_view kFirstBindingToken{ "GAUS" };
    static constexpr std::array<EnumMappingExpectation<rg::PainterType>, 3> kExpectations{{
        { "gaus", rg::PainterType::GAUS },
        { "0", rg::PainterType::GAUS },
        { "atom", rg::PainterType::ATOM },
    }};
};

template <>
struct EnumMappingTraits<rg::PrinterType>
{
    static constexpr std::string_view kFirstBindingToken{ "ATOM_POSITION" };
    static constexpr std::array<EnumMappingExpectation<rg::PrinterType>, 2> kExpectations{{
        { "atom_out", rg::PrinterType::ATOM_OUTLIER },
        { "3", rg::PrinterType::ATOM_OUTLIER },
    }};
};

template <>
struct EnumMappingTraits<rg::PotentialModel>
{
    static constexpr std::string_view kFirstBindingToken{ "SINGLE_GAUS" };
    static constexpr std::array<EnumMappingExpectation<rg::PotentialModel>, 2> kExpectations{{
        { "five", rg::PotentialModel::FIVE_GAUS_CHARGE },
        { "1", rg::PotentialModel::FIVE_GAUS_CHARGE },
    }};
};

template <>
struct EnumMappingTraits<rg::PartialCharge>
{
    static constexpr std::string_view kFirstBindingToken{ "NEUTRAL" };
    static constexpr std::array<EnumMappingExpectation<rg::PartialCharge>, 2> kExpectations{{
        { "amber", rg::PartialCharge::AMBER },
        { "2", rg::PartialCharge::AMBER },
    }};
};

template <>
struct EnumMappingTraits<rg::TesterType>
{
    static constexpr std::string_view kFirstBindingToken{ "BENCHMARK" };
    static constexpr std::array<EnumMappingExpectation<rg::TesterType>, 2> kExpectations{{
        { "benchmark", rg::TesterType::BENCHMARK },
        { "0", rg::TesterType::BENCHMARK },
    }};
};

template <typename EnumType>
void AssertEnumMappingsStayInSync()
{
    const auto cli_map{ rg::internal::BuildCommandEnumCliMap<EnumType>() };
    for (const auto & expectation : EnumMappingTraits<EnumType>::kExpectations)
    {
        EXPECT_EQ(cli_map.at(std::string(expectation.token)), expectation.value);
    }

    const auto binding_entries{ rg::internal::GetCommandEnumBindingEntries<EnumType>() };
    ASSERT_FALSE(binding_entries.empty());
    EXPECT_EQ(binding_entries.front().token, EnumMappingTraits<EnumType>::kFirstBindingToken);
    EXPECT_EQ(BuildValueSetFromCLI<EnumType>(), BuildValueSetFromBindings<EnumType>());
}

template <typename EnumType>
class CommandEnumMappingTest : public testing::Test {};

using CommandEnumTypes = testing::Types<
    rg::PainterType,
    rg::PrinterType,
    rg::PotentialModel,
    rg::PartialCharge,
    rg::TesterType>;

TYPED_TEST_SUITE(CommandEnumMappingTest, CommandEnumTypes, );

} // namespace

TEST(PublicHeaderSurfaceTest, CorePublicHeadersMatchApprovedSurface) {
    const std::vector<std::string> expected{
        "core/command/CommandApi.hpp",
        "core/command/CommandEnums.hpp",
        "core/command/CommandManifest.def",
        "core/painter/AtomPainter.hpp",
        "core/painter/ComparisonPainter.hpp",
        "core/painter/DemoPainter.hpp",
        "core/painter/GausPainter.hpp",
        "core/painter/ModelPainter.hpp",
        "core/painter/PainterBase.hpp"};

    EXPECT_EQ(contract_test_support::CollectPublicHeadersForDomain("core"), expected);
}

TEST(PublicHeaderSurfaceTest, CommandApiExposesStableDefaultDatabasePathHelper) {
    const auto default_path{ rhbm_gem::GetDefaultDatabasePath() };
    EXPECT_FALSE(default_path.empty());
    EXPECT_EQ(default_path.filename(), "database.sqlite");
}

TYPED_TEST(CommandEnumMappingTest, CliAndBindingMappingsStayInSync)
{
    AssertEnumMappingsStayInSync<TypeParam>();
}
