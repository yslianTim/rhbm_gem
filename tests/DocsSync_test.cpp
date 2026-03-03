#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "BuiltInCommandCatalog.hpp"

namespace rg = rhbm_gem;

namespace {

constexpr std::string_view kManifestBeginMarker{
    "<!-- BEGIN GENERATED: built-in-command-manifest -->"
};
constexpr std::string_view kManifestEndMarker{
    "<!-- END GENERATED: built-in-command-manifest -->"
};
constexpr std::string_view kMatrixBeginMarker{
    "<!-- BEGIN GENERATED: command-surface-matrix -->"
};
constexpr std::string_view kMatrixEndMarker{
    "<!-- END GENERATED: command-surface-matrix -->"
};
constexpr std::string_view kPythonBeginMarker{
    "<!-- BEGIN GENERATED: python-public-command-surface -->"
};
constexpr std::string_view kPythonEndMarker{
    "<!-- END GENERATED: python-public-command-surface -->"
};

std::string YesNo(bool value)
{
    return value ? "yes" : "no";
}

std::string RenderManifest(const std::vector<rg::CommandDescriptor> & commands)
{
    std::ostringstream output;
    std::size_t index{ 1 };
    for (const auto & info : commands)
    {
        output << index++ << ". `" << info.name << "`\n";
    }
    return output.str();
}

std::string RenderSurfaceMatrix(const std::vector<rg::CommandDescriptor> & commands)
{
    std::ostringstream output;
    output
        << "| Command | Uses database at runtime | Uses output folder | Exposed to Python | Hidden deprecated `--database` alias |\n"
        << "| --- | --- | --- | --- | --- |\n";
    for (const auto & info : commands)
    {
        output
            << "| `" << info.name << "`"
            << " | " << YesNo(rg::UsesDatabaseAtRuntime(info.database_usage))
            << " | " << YesNo(rg::UsesOutputFolder(info.surface))
            << " | " << YesNo(rg::IsPythonPublic(info.binding_exposure))
            << " | " << YesNo(rg::HasCommonOption(
                info.surface.deprecated_hidden_options,
                rg::CommonOption::Database))
            << " |\n";
    }
    return output.str();
}

std::string RenderPythonPublicSurface(const std::vector<rg::CommandDescriptor> & commands)
{
    static constexpr std::array<std::string_view, 3> kDiagnosticsTypes{
        "LogLevel",
        "ValidationPhase",
        "ValidationIssue"
    };
    static constexpr std::array<std::string_view, 3> kDiagnosticsMethods{
        "PrepareForExecution()",
        "HasValidationErrors()",
        "GetValidationIssues()"
    };

    std::ostringstream output;
    output << "### Python-public command classes\n";
    for (const auto & descriptor : commands)
    {
        if (!rg::IsPythonPublic(descriptor.binding_exposure)) continue;
        output << "- `" << descriptor.python_binding_name << "`\n";
    }

    output << "\n### Shared diagnostics types\n";
    for (const auto type_name : kDiagnosticsTypes)
    {
        output << "- `" << type_name << "`\n";
    }

    output << "\n### Shared diagnostics methods on Python-public commands\n";
    for (const auto method_name : kDiagnosticsMethods)
    {
        output << "- `" << method_name << "`\n";
    }

    return output.str();
}

std::string ReadGeneratedBlock(
    const std::string & content,
    std::string_view begin_marker,
    std::string_view end_marker)
{
    const auto begin_pos{ content.find(begin_marker) };
    EXPECT_NE(begin_pos, std::string::npos);
    const auto end_pos{ content.find(end_marker) };
    EXPECT_NE(end_pos, std::string::npos);
    EXPECT_LT(begin_pos, end_pos);
    if (begin_pos == std::string::npos || end_pos == std::string::npos || begin_pos >= end_pos)
    {
        return {};
    }
    const auto content_start{ begin_pos + begin_marker.size() };
    return content.substr(content_start, end_pos - content_start);
}

} // namespace

TEST(DocsSyncTest, GeneratedBlocksMatchBuiltInCommandCatalog)
{
    const auto & commands{ rg::BuiltInCommandCatalog() };
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

    EXPECT_EQ(
        ReadGeneratedBlock(doc_content, kManifestBeginMarker, kManifestEndMarker),
        "\n" + RenderManifest(commands));
    EXPECT_EQ(
        ReadGeneratedBlock(doc_content, kMatrixBeginMarker, kMatrixEndMarker),
        "\n" + RenderSurfaceMatrix(commands));
    EXPECT_EQ(
        ReadGeneratedBlock(doc_content, kPythonBeginMarker, kPythonEndMarker),
        "\n" + RenderPythonPublicSurface(commands));

    EXPECT_NE(
        doc_content.find("Application callbacks invoke only `Execute()`."),
        std::string::npos);
    EXPECT_NE(
        doc_content.find(
            "`Execute()` internally decides whether `PrepareForExecution()` must run."),
        std::string::npos);
}
