#include <gtest/gtest.h>

#include <filesystem>
#include <string>
#include <CLI/CLI.hpp>

#include "internal/CommandOptionBinding.hpp"
#include <rhbm_gem/core/command/OptionEnumClass.hpp>

namespace rg = rhbm_gem;

TEST(CommandCLIHelperTest, AddScalarOptionShowsDefaultValueInHelpAndParsesInput) {
    CLI::App app{"helper-test"};
    int count{0};
    rg::command_cli::AddScalarOption<int>(
        &app,
        "--count",
        [&](int value) { count = value; },
        "Count option",
        42);

    const std::string help_text{app.help()};
    EXPECT_NE(help_text.find("--count"), std::string::npos);
    EXPECT_NE(help_text.find("42"), std::string::npos);

    app.parse("--count 7", false);
    EXPECT_EQ(count, 7);
}

TEST(CommandCLIHelperTest, AddPathOptionPreservesRequiredBehavior) {
    CLI::App app{"helper-test"};
    std::filesystem::path parsed_path;
    rg::command_cli::AddPathOption(
        &app,
        "--input",
        [&](const std::filesystem::path& value) { parsed_path = value; },
        "Input path",
        std::nullopt,
        true);

    EXPECT_THROW(app.parse("", false), CLI::RequiredError);

    CLI::App parse_app{"helper-test"};
    rg::command_cli::AddPathOption(
        &parse_app,
        "--input",
        [&](const std::filesystem::path& value) { parsed_path = value; },
        "Input path",
        std::nullopt,
        true);
    parse_app.parse("--input fixture.map", false);
    EXPECT_EQ(parsed_path, std::filesystem::path("fixture.map"));
}

TEST(CommandCLIHelperTest, AddEnumOptionUsesEnumTransformer) {
    CLI::App app{"helper-test"};
    rg::PrinterType printer{rg::PrinterType::ATOM_POSITION};
    rg::command_cli::AddEnumOption<rg::PrinterType>(
        &app,
        "--printer",
        [&](rg::PrinterType value) { printer = value; },
        "Printer option",
        std::nullopt,
        true);

    app.parse("--printer map", false);
    EXPECT_EQ(printer, rg::PrinterType::MAP_VALUE);
}
