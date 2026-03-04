#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>

#include "BuiltInCommandCatalogInternal.hpp"

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
    "<!-- BEGIN GENERATED: built-in-python-command-surface -->"
};
constexpr std::string_view kPythonEndMarker{
    "<!-- END GENERATED: built-in-python-command-surface -->"
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
        << "| Command | Uses database at runtime | Uses output folder |\n"
        << "| --- | --- | --- |\n";
    for (const auto & info : commands)
    {
        output
            << "| `" << info.name << "`"
            << " | " << YesNo(rg::UsesDatabaseAtRuntime(info.common_options))
            << " | " << YesNo(rg::UsesOutputFolder(info.common_options))
            << " |\n";
    }
    return output.str();
}

std::string RenderBuiltInPythonSurface(const std::vector<rg::CommandDescriptor> & commands)
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
    output << "### Built-in Python command classes\n";
    for (const auto & descriptor : commands)
    {
        output << "- `" << descriptor.python_binding_name << "`\n";
    }

    output << "\n### Shared diagnostics types\n";
    for (const auto type_name : kDiagnosticsTypes)
    {
        output << "- `" << type_name << "`\n";
    }

    output << "\n### Shared diagnostics methods on built-in Python commands\n";
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
        "\n" + RenderBuiltInPythonSurface(commands));

    EXPECT_NE(
        doc_content.find("Application callbacks invoke only `Execute()`."),
        std::string::npos);
    EXPECT_NE(
        doc_content.find(
            "`Execute()` internally decides whether `PrepareForExecution()` must run."),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("CommandRegistrar<"),
        std::string::npos);
    EXPECT_NE(
        doc_content.find(
            "The project does not currently provide a self-registration API for commands."),
        std::string::npos);
    EXPECT_NE(
        doc_content.find(
            "Database usage is derived directly from whether the built-in descriptor's `common_options` mask includes `CommonOption::Database`."),
        std::string::npos);
    EXPECT_NE(
        doc_content.find(
            "All built-in commands must provide a Python binding name in the built-in catalog."),
        std::string::npos);
    EXPECT_NE(
        doc_content.find(
            "The built-in catalog is an internal registration mechanism, not an installed public C++ API."),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("DatabaseUsage"),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("BindingExposure"),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("Python-public command classes"),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("subset of commands through `bindings/CoreBindings.cpp`"),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("Exposed to Python"),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("include/core/BuiltInCommandCatalog.hpp"),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("CommandSurface"),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("ProcessTypedFile"),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("OptionalProcessTypedFile"),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find("LoadTypedObject"),
        std::string::npos);
    EXPECT_NE(
        doc_content.find("src/core/BuiltInCommandCatalogInternal.hpp"),
        std::string::npos);
}
