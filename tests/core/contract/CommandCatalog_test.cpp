#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "command/detail/CommandCatalog.hpp"
#include "command/detail/CommandEnumCatalog.hpp"
#include <rhbm_gem/core/command/CommandSystem.hpp>

using namespace rhbm_gem;

namespace {

template <typename EnumType>
std::map<std::string, EnumType> BuildCliMapFromTraits()
{
    std::map<std::string, EnumType> option_map;
    for (const auto & option : command_internal::CommandEnumTraits<EnumType>::kOptions)
    {
        for (const auto alias : option.cli_aliases)
        {
            option_map.emplace(std::string(alias), option.value);
        }
    }
    return option_map;
}

template <typename EnumType>
std::set<int> BuildCliValueSetFromTraits()
{
    std::set<int> values;
    for (const auto & option : command_internal::CommandEnumTraits<EnumType>::kOptions)
    {
        for (const auto alias : option.cli_aliases)
        {
            static_cast<void>(alias);
            values.insert(static_cast<int>(option.value));
        }
    }
    return values;
}

template <typename EnumType>
std::set<int> BuildBindingValueSetFromTraits()
{
    std::set<int> values;
    for (const auto & option : command_internal::CommandEnumTraits<EnumType>::kOptions)
    {
        values.insert(static_cast<int>(option.value));
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
struct EnumMappingTraits<PainterType>
{
    static constexpr std::string_view kFirstBindingToken{ "GAUS" };
    static constexpr std::array<EnumMappingExpectation<PainterType>, 3> kExpectations{{
        { "gaus", PainterType::GAUS },
        { "0", PainterType::GAUS },
        { "atom", PainterType::ATOM },
    }};
};

template <>
struct EnumMappingTraits<PrinterType>
{
    static constexpr std::string_view kFirstBindingToken{ "ATOM_POSITION" };
    static constexpr std::array<EnumMappingExpectation<PrinterType>, 2> kExpectations{{
        { "atom_out", PrinterType::ATOM_OUTLIER },
        { "3", PrinterType::ATOM_OUTLIER },
    }};
};

template <>
struct EnumMappingTraits<PotentialModel>
{
    static constexpr std::string_view kFirstBindingToken{ "SINGLE_GAUS" };
    static constexpr std::array<EnumMappingExpectation<PotentialModel>, 2> kExpectations{{
        { "five", PotentialModel::FIVE_GAUS_CHARGE },
        { "1", PotentialModel::FIVE_GAUS_CHARGE },
    }};
};

template <>
struct EnumMappingTraits<PartialCharge>
{
    static constexpr std::string_view kFirstBindingToken{ "NEUTRAL" };
    static constexpr std::array<EnumMappingExpectation<PartialCharge>, 2> kExpectations{{
        { "amber", PartialCharge::AMBER },
        { "2", PartialCharge::AMBER },
    }};
};

template <>
struct EnumMappingTraits<TesterType>
{
    static constexpr std::string_view kFirstBindingToken{ "BENCHMARK" };
    static constexpr std::array<EnumMappingExpectation<TesterType>, 4> kExpectations{{
        { "benchmark", TesterType::BENCHMARK },
        { "0", TesterType::BENCHMARK },
        { "neighbor_distance", TesterType::NEIGHBOR_DISTANCE },
        { "5", TesterType::NEIGHBOR_DISTANCE },
    }};
};

template <typename EnumType>
void AssertEnumMappingsStayInSync()
{
    const auto cli_map{ BuildCliMapFromTraits<EnumType>() };
    for (const auto & expectation : EnumMappingTraits<EnumType>::kExpectations)
    {
        EXPECT_EQ(cli_map.at(std::string(expectation.token)), expectation.value);
    }

    const auto & enum_options{ command_internal::CommandEnumTraits<EnumType>::kOptions };
    ASSERT_FALSE(enum_options.empty());
    EXPECT_EQ(enum_options.front().binding_token, EnumMappingTraits<EnumType>::kFirstBindingToken);
    EXPECT_EQ(BuildCliValueSetFromTraits<EnumType>(), BuildBindingValueSetFromTraits<EnumType>());
}

template <typename EnumType>
class CommandEnumMappingTest : public testing::Test {};

using CommandEnumTypes = testing::Types<
    PainterType,
    PrinterType,
    PotentialModel,
    PartialCharge,
    TesterType>;

TYPED_TEST_SUITE(CommandEnumMappingTest, CommandEnumTypes, );

std::vector<std::string> BuildExpectedCommandNames()
{
    std::vector<std::string> expected{
        "potential_analysis",
        "potential_display",
        "result_dump",
        "map_simulation",
        "rhbm_test",
    };
#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    expected.emplace_back("map_visualization");
    expected.emplace_back("position_estimation");
#endif
    return expected;
}

std::vector<std::string> BuildCommandNames()
{
    std::vector<std::string> names;
    command_internal::VisitCommandCatalog([&](const auto & entry)
    {
        names.emplace_back(entry.cli_name);
    });
    return names;
}

} // namespace

TEST(CommandCatalogTest, CommandCatalogKeepsExpectedOrder)
{
    EXPECT_EQ(BuildExpectedCommandNames(), BuildCommandNames());
}

TYPED_TEST(CommandEnumMappingTest, CliAndBindingMappingsStayInSync)
{
    AssertEnumMappingsStayInSync<TypeParam>();
}

TEST(CommandCatalogTest, ExperimentalCommandVisibilityFollowsFeatureGuard)
{
    const std::vector<std::string> expected_experimental_names{
        "map_visualization",
        "position_estimation",
    };
    const auto command_names{ BuildCommandNames() };

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    ASSERT_GE(command_names.size(), expected_experimental_names.size());
    const std::vector<std::string> trailing_names{
        command_names.end() - static_cast<std::ptrdiff_t>(expected_experimental_names.size()),
        command_names.end()
    };
    EXPECT_EQ(trailing_names, expected_experimental_names);
#else
    for (const auto & experimental_name : expected_experimental_names)
    {
        EXPECT_EQ(
            std::find(command_names.begin(), command_names.end(), experimental_name),
            command_names.end());
    }
#endif
}

TEST(CommandCatalogTest, RunCommandCLIReturnsSuccessForHelpRequest)
{
    char program[]{ "RHBM-GEM" };
    char help_flag[]{ "--help" };
    char * argv[]{ program, help_flag };
    const int argc{ static_cast<int>(std::size(argv)) };

    EXPECT_EQ(RunCommandCLI(argc, argv), 0);
}

TEST(CommandCatalogTest, RunCommandCLIReturnsFailureForMissingCommandInput)
{
    char program[]{ "RHBM-GEM" };
    char command[]{ "map_simulation" };
    char * argv[]{ program, command };
    const int argc{ static_cast<int>(std::size(argv)) };

    EXPECT_NE(RunCommandCLI(argc, argv), 0);
}
