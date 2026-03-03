#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "CommandRegistry.hpp"
#include "RegisterBuiltInCommands.hpp"

namespace rg = rhbm_gem;

namespace {

constexpr std::string_view kBeginMarker{
    "<!-- BEGIN GENERATED: command-surface-matrix -->"
};
constexpr std::string_view kEndMarker{
    "<!-- END GENERATED: command-surface-matrix -->"
};

std::string YesNo(bool value)
{
    return value ? "yes" : "no";
}

std::string RenderSurfaceMatrix(const std::vector<rg::CommandRegistry::CommandInfo> & commands)
{
    std::ostringstream output;
    output
        << "| Command | Uses database at runtime | Uses output folder | Exposed to Python | Hidden deprecated `--database` alias |\n"
        << "| --- | --- | --- | --- | --- |\n";
    for (const auto & info : commands)
    {
        output
            << "| `" << info.name << "`"
            << " | " << YesNo(info.surface.requires_database_manager)
            << " | " << YesNo(info.surface.uses_output_folder)
            << " | " << YesNo(info.surface.exposed_to_python)
            << " | " << YesNo(HasCommonOption(
                info.surface.deprecated_hidden_options,
                rg::CommonOption::Database))
            << " |\n";
    }
    return output.str();
}

} // namespace

TEST(DocsSyncTest, CommandSurfaceMatrixMatchesRegistryMetadata)
{
    rg::RegisterBuiltInCommands();
    const auto & commands{ rg::CommandRegistry::Instance().GetCommands() };
    const auto doc_path{
        std::filesystem::path(__FILE__).parent_path().parent_path()
        / "docs" / "developer" / "architecture" / "command-architecture.md"
    };

    std::ifstream doc_stream(doc_path);
    ASSERT_TRUE(doc_stream.is_open());
    const std::string doc_content{
        std::istreambuf_iterator<char>(doc_stream),
        std::istreambuf_iterator<char>()
    };

    const auto begin_pos{ doc_content.find(kBeginMarker) };
    ASSERT_NE(begin_pos, std::string::npos);
    const auto end_pos{ doc_content.find(kEndMarker) };
    ASSERT_NE(end_pos, std::string::npos);
    ASSERT_LT(begin_pos, end_pos);

    const auto content_start{ begin_pos + kBeginMarker.size() };
    const std::string matrix_block{
        doc_content.substr(content_start, end_pos - content_start)
    };

    const std::string expected_block{ "\n" + RenderSurfaceMatrix(commands) };
    EXPECT_EQ(matrix_block, expected_block);
}
