#include <gtest/gtest.h>

#include <rhbm_gem/core/CommandSystem.hpp>

#include "support/CommandTestHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace rhbm_gem::core;

namespace {

constexpr std::string_view kInputHeader{
    "serial id,residue,spot,neighbor count,"
    "signal peeling ratio,tail peeling ratio,"
    "amplitude 1st,amplitude 2nd,amplitude 3rd,"
    "width 1st,width 2nd,width 3rd,"
    "offset 1st,offset 2nd,offset 3rd,"
    "amplitude rank 1st,amplitude rank 2nd,amplitude rank 3rd,"
    "width rank 1st,width rank 2nd,width rank 3rd,"
    "offset rank 1st,offset rank 2nd,offset rank 3rd"
};
constexpr std::string_view kOldInputHeader{
    "serial id,residue,spot,neighbor count,peeling ratio,"
    "amplitude 1st,amplitude 2nd,amplitude 3rd,"
    "width 1st,width 2nd,width 3rd,"
    "offset 1st,offset 2nd,offset 3rd,"
    "amplitude rank 1st,amplitude rank 2nd,amplitude rank 3rd,"
    "width rank 1st,width rank 2nd,width rank 3rd,"
    "offset rank 1st,offset rank 2nd,offset rank 3rd"
};
constexpr std::size_t kFeatureCount{ 21 };
constexpr std::size_t kInputColumnCount{ 24 };
constexpr std::size_t kOutputColumnCount{ 26 };

std::vector<std::string> BuildFields(
    std::size_t observation,
    std::optional<std::size_t> constant_feature = std::nullopt,
    bool apply_feature_specific_affine_transform = false)
{
    std::vector<std::string> fields{
        std::to_string(observation + 1),
        "ALA",
        "CA",
    };
    fields.reserve(kInputColumnCount);
    for (std::size_t feature = 0; feature < kFeatureCount; ++feature)
    {
        const auto left{ (observation + 1) * (feature + 2) };
        const auto right{ (observation * observation + 3 * feature) % (feature + 3) };
        double value{ static_cast<double>(left + right) };
        if (constant_feature && feature == *constant_feature)
        {
            value = 7;
        }
        if (apply_feature_specific_affine_transform)
        {
            value = value * static_cast<double>(feature + 2)
                + static_cast<double>(17 * feature);
        }
        std::ostringstream text;
        text << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
        fields.push_back(text.str());
    }
    return fields;
}

std::string JoinFields(const std::vector<std::string> & fields)
{
    std::string row;
    for (std::size_t index = 0; index < fields.size(); ++index)
    {
        if (index != 0) row += ',';
        row += fields[index];
    }
    return row;
}

std::vector<std::string> MakeRows(
    std::size_t count,
    std::optional<std::size_t> constant_feature = std::nullopt,
    bool apply_feature_specific_affine_transform = false)
{
    std::vector<std::string> rows;
    rows.reserve(count);
    for (std::size_t observation = 0; observation < count; ++observation)
    {
        rows.push_back(JoinFields(BuildFields(
            observation,
            constant_feature,
            apply_feature_specific_affine_transform)));
    }
    return rows;
}

void WriteCsv(
    const std::filesystem::path & path,
    const std::vector<std::string> & rows,
    std::string_view header = kInputHeader,
    std::string_view line_ending = "\n")
{
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output{ path, std::ios::binary };
    ASSERT_TRUE(output.is_open());
    output << header << line_ending;
    for (const auto & row : rows)
    {
        output << row << line_ending;
    }
    ASSERT_TRUE(output.good());
}

std::vector<std::string> ReadLines(const std::filesystem::path & path)
{
    std::ifstream input{ path };
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line))
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }
    return lines;
}

std::vector<std::string> SplitFields(std::string_view row)
{
    std::vector<std::string> fields;
    std::size_t start{ 0 };
    while (true)
    {
        const auto delimiter{ row.find(',', start) };
        if (delimiter == std::string_view::npos)
        {
            fields.emplace_back(row.substr(start));
            return fields;
        }
        fields.emplace_back(row.substr(start, delimiter - start));
        start = delimiter + 1;
    }
}

bool HasIssue(
    const CommandResult & result,
    std::string_view option_name,
    std::string_view message_fragment = {})
{
    return std::any_of(
        result.issues.begin(),
        result.issues.end(),
        [option_name, message_fragment](const CommandDiagnostic & issue)
        {
            return issue.option_name == option_name
                && (message_fragment.empty()
                    || issue.message.find(message_fragment) != std::string::npos);
        });
}

UmapEmbeddingRequest MakeRequest(
    const std::filesystem::path & input_path,
    const std::filesystem::path & output_dir)
{
    UmapEmbeddingRequest request;
    request.input_csv_path = input_path;
    request.output_dir = output_dir;
    request.num_neighbors = 4;
    request.num_epochs = 10;
    request.random_seed = 1234;
    request.job_count = 1;
    return request;
}

std::vector<std::array<double, 2>> ReadCoordinates(
    const std::filesystem::path & output_path)
{
    const auto lines{ ReadLines(output_path) };
    std::vector<std::array<double, 2>> coordinates;
    for (std::size_t line = 1; line < lines.size(); ++line)
    {
        const auto fields{ SplitFields(lines[line]) };
        if (fields.size() != kOutputColumnCount) return {};
        coordinates.push_back({
            std::stod(fields[kInputColumnCount]),
            std::stod(fields[kInputColumnCount + 1])
        });
    }
    return coordinates;
}

} // namespace

TEST(UmapEmbeddingCommandTest, DefaultsMatchThePublicContract)
{
    UmapEmbeddingRequest request;

    EXPECT_EQ(request.num_neighbors, 15);
    EXPECT_DOUBLE_EQ(request.min_dist, 0.1);
    EXPECT_EQ(request.num_epochs, 0);
    EXPECT_EQ(request.random_seed, 42u);
}

TEST(UmapEmbeddingCommandTest, ProducesTwoFiniteCoordinatesAndPreservesCrLfRows)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_success" };
    const auto input_path{ temp_dir.path() / "local_fitting_result_case.csv" };
    const auto output_dir{ temp_dir.path() / "output" };
    const auto input_rows{ MakeRows(8, 13) };
    WriteCsv(input_path, input_rows, kInputHeader, "\r\n");

    auto request{ MakeRequest(input_path, output_dir) };
    request.num_neighbors = 15;
    const auto result{ RunCommand(request) };

    ASSERT_TRUE(result.succeeded);
    EXPECT_TRUE(HasIssue(result, "-i,--input", "amplitude rank 2nd"));
    EXPECT_TRUE(HasIssue(result, "--neighbors", "limited to 7"));

    const auto output_path{ output_dir / "umap_embedding_case.csv" };
    const auto output_lines{ ReadLines(output_path) };
    ASSERT_EQ(output_lines.size(), input_rows.size() + 1);
    EXPECT_EQ(output_lines.front(), std::string(kInputHeader) + ",umap x,umap y");
    for (std::size_t row = 0; row < input_rows.size(); ++row)
    {
        const auto fields{ SplitFields(output_lines[row + 1]) };
        ASSERT_EQ(fields.size(), kOutputColumnCount);
        EXPECT_EQ(
            output_lines[row + 1].substr(0, input_rows[row].size()),
            input_rows[row]);
        EXPECT_TRUE(std::isfinite(std::stod(fields[kInputColumnCount])));
        EXPECT_TRUE(std::isfinite(std::stod(fields[kInputColumnCount + 1])));
    }
}

TEST(UmapEmbeddingCommandTest, RejectsInvalidOrReorderedHeaderAndDoesNotCreateOutput)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_header" };
    const std::vector<std::string> invalid_headers{
        "serial id,residue,spot",
        std::string{ kOldInputHeader },
        "serial id,spot,residue,neighbor count,"
            "signal peeling ratio,tail peeling ratio,"
            "amplitude 1st,amplitude 2nd,amplitude 3rd,"
            "width 1st,width 2nd,width 3rd,"
            "offset 1st,offset 2nd,offset 3rd,"
            "amplitude rank 1st,amplitude rank 2nd,amplitude rank 3rd,"
            "width rank 1st,width rank 2nd,width rank 3rd,"
            "offset rank 1st,offset rank 2nd,offset rank 3rd",
    };
    for (std::size_t index = 0; index < invalid_headers.size(); ++index)
    {
        SCOPED_TRACE(index);
        const auto case_dir{ temp_dir.path() / std::to_string(index) };
        const auto input_path{ case_dir / "local_fitting_result_bad.csv" };
        const auto output_dir{ case_dir / "output" };
        WriteCsv(input_path, MakeRows(4), invalid_headers[index]);

        const auto result{ RunCommand(MakeRequest(input_path, output_dir)) };

        EXPECT_FALSE(result.succeeded);
        EXPECT_TRUE(HasIssue(result, "-i,--input", "line 1, column 'header'"));
        EXPECT_FALSE(std::filesystem::exists(output_dir / "umap_embedding_bad.csv"));
    }
}

TEST(UmapEmbeddingCommandTest, RejectsMalformedAndNonFiniteFeatureValues)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_malformed" };
    struct InvalidCase
    {
        std::string name;
        std::vector<std::string> fields;
        std::string diagnostic;
    };

    auto missing_field{ BuildFields(0) };
    missing_field.pop_back();
    auto empty_field{ BuildFields(0) };
    empty_field[3] = "";
    auto invalid_signal_ratio{ BuildFields(0) };
    invalid_signal_ratio[4] = "12x";
    auto nan_tail_ratio{ BuildFields(0) };
    nan_tail_ratio[5] = "nan";
    auto infinity_value{ BuildFields(0) };
    infinity_value[6] = "inf";
    auto rank_nan_value{ BuildFields(0) };
    rank_nan_value[15] = "nan";

    const std::vector<InvalidCase> cases{
        { "missing", std::move(missing_field), "column 'layout'" },
        { "empty", std::move(empty_field), "column 'neighbor count'" },
        { "invalid_signal", std::move(invalid_signal_ratio), "column 'signal peeling ratio'" },
        { "nan_tail", std::move(nan_tail_ratio), "column 'tail peeling ratio'" },
        { "infinity", std::move(infinity_value), "column 'amplitude 1st'" },
        { "rank_nan", std::move(rank_nan_value), "column 'amplitude rank 1st'" },
    };

    for (const auto & invalid_case : cases)
    {
        SCOPED_TRACE(invalid_case.name);
        const auto case_dir{ temp_dir.path() / invalid_case.name };
        const auto input_path{ case_dir / "local_fitting_result_case.csv" };
        auto rows{ MakeRows(4) };
        rows.front() = JoinFields(invalid_case.fields);
        WriteCsv(input_path, rows);

        const auto result{ RunCommand(MakeRequest(input_path, case_dir / "output")) };

        EXPECT_FALSE(result.succeeded);
        EXPECT_TRUE(HasIssue(result, "-i,--input", "line 2"));
        EXPECT_TRUE(HasIssue(result, "-i,--input", invalid_case.diagnostic));
        EXPECT_FALSE(std::filesystem::exists(case_dir / "output" / "umap_embedding_case.csv"));
    }
}

TEST(UmapEmbeddingCommandTest, RejectsQuotedRowsAndInputsWithFewerThanThreeRows)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_row_contract" };

    const auto quoted_path{ temp_dir.path() / "quoted.csv" };
    auto quoted_rows{ MakeRows(3) };
    quoted_rows.front().replace(quoted_rows.front().find("ALA"), 3, "\"ALA\"");
    WriteCsv(quoted_path, quoted_rows);
    const auto quoted_result{
        RunCommand(MakeRequest(quoted_path, temp_dir.path() / "quoted_output"))
    };
    EXPECT_FALSE(quoted_result.succeeded);
    EXPECT_TRUE(HasIssue(quoted_result, "-i,--input", "quoted CSV fields"));

    const auto short_path{ temp_dir.path() / "short.csv" };
    WriteCsv(short_path, MakeRows(2));
    const auto short_result{
        RunCommand(MakeRequest(short_path, temp_dir.path() / "short_output"))
    };
    EXPECT_FALSE(short_result.succeeded);
    EXPECT_TRUE(HasIssue(short_result, "-i,--input", "at least 3 data rows"));
}

TEST(UmapEmbeddingCommandTest, RejectsAllConstantFeaturesWithoutWritingOutput)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_all_constant" };
    const auto input_path{ temp_dir.path() / "local_fitting_result_constant.csv" };
    const auto output_dir{ temp_dir.path() / "output" };
    std::vector<std::string> rows;
    for (std::size_t observation = 0; observation < 4; ++observation)
    {
        rows.push_back(JoinFields(BuildFields(0)));
    }
    WriteCsv(input_path, rows);

    const auto result{ RunCommand(MakeRequest(input_path, output_dir)) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasIssue(
        result,
        "-i,--input",
        "All 11 selected UMAP feature columns are constant"));
    EXPECT_FALSE(std::filesystem::exists(output_dir / "umap_embedding_constant.csv"));
}

TEST(UmapEmbeddingCommandTest, UsesRankColumnsAsEmbeddingFeatures)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_rank_features" };
    const auto input_path{ temp_dir.path() / "local_fitting_result_rank.csv" };
    const auto output_dir{ temp_dir.path() / "output" };
    std::vector<std::string> rows;
    for (std::size_t observation = 0; observation < 6; ++observation)
    {
        auto fields{ BuildFields(0) };
        fields[0] = std::to_string(observation + 1);
        for (const auto feature : { 13u, 16u, 19u })
        {
            fields[3 + feature] = std::to_string((observation + feature) % 4 + 1);
        }
        rows.push_back(JoinFields(fields));
    }
    WriteCsv(input_path, rows);

    const auto result{ RunCommand(MakeRequest(input_path, output_dir)) };

    EXPECT_TRUE(result.succeeded);
    EXPECT_TRUE(std::filesystem::exists(output_dir / "umap_embedding_rank.csv"));
}

TEST(UmapEmbeddingCommandTest, UsesBothPeelingRatioColumnsAsEmbeddingFeatures)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_peeling_ratio_features" };
    for (const auto ratio_feature : { 1u, 2u })
    {
        SCOPED_TRACE(ratio_feature);
        const auto case_dir{ temp_dir.path() / std::to_string(ratio_feature) };
        const auto input_path{ case_dir / "local_fitting_result_ratio.csv" };
        const auto output_dir{ case_dir / "output" };
        std::vector<std::string> rows;
        for (std::size_t observation = 0; observation < 6; ++observation)
        {
            auto fields{ BuildFields(0) };
            fields[0] = std::to_string(observation + 1);
            fields[3 + ratio_feature] = std::to_string(observation + 1);
            rows.push_back(JoinFields(fields));
        }
        WriteCsv(input_path, rows);

        const auto result{ RunCommand(MakeRequest(input_path, output_dir)) };

        EXPECT_TRUE(result.succeeded);
        EXPECT_TRUE(std::filesystem::exists(output_dir / "umap_embedding_ratio.csv"));
    }
}

TEST(UmapEmbeddingCommandTest, RejectsInvalidOptionsAndNonRegularInput)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_options" };
    const auto input_path{ temp_dir.path() / "input.csv" };
    WriteCsv(input_path, MakeRows(4));

    auto neighbors_request{ MakeRequest(input_path, temp_dir.path() / "neighbors") };
    neighbors_request.num_neighbors = 1;
    EXPECT_TRUE(HasIssue(RunCommand(neighbors_request), "--neighbors"));

    auto min_dist_request{ MakeRequest(input_path, temp_dir.path() / "min_dist") };
    min_dist_request.min_dist = 1.01;
    EXPECT_TRUE(HasIssue(RunCommand(min_dist_request), "--min-dist"));

    min_dist_request.min_dist = -0.01;
    EXPECT_TRUE(HasIssue(RunCommand(min_dist_request), "--min-dist"));

    min_dist_request.min_dist = std::numeric_limits<double>::quiet_NaN();
    EXPECT_TRUE(HasIssue(RunCommand(min_dist_request), "--min-dist"));

    auto epochs_request{ MakeRequest(input_path, temp_dir.path() / "epochs") };
    epochs_request.num_epochs = -1;
    EXPECT_TRUE(HasIssue(RunCommand(epochs_request), "--epochs"));

    auto directory_request{ MakeRequest(temp_dir.path(), temp_dir.path() / "directory") };
    const auto directory_result{ RunCommand(directory_request) };
    EXPECT_FALSE(directory_result.succeeded);
    EXPECT_TRUE(HasIssue(directory_result, "-i,--input", "regular file"));
}

TEST(UmapEmbeddingCommandTest, AcceptsDocumentedOptionBoundaries)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_boundaries" };
    const auto input_path{ temp_dir.path() / "custom.csv" };
    WriteCsv(input_path, MakeRows(5));
    auto request{ MakeRequest(input_path, temp_dir.path() / "output") };
    request.num_neighbors = 2;
    request.min_dist = 0;
    request.num_epochs = 1;

    ASSERT_TRUE(RunCommand(request).succeeded);
    EXPECT_TRUE(std::filesystem::exists(
        request.output_dir / "umap_embedding_custom.csv"));

    request.output_dir = temp_dir.path() / "output_max";
    request.min_dist = 1;
    request.num_epochs = 0;
    EXPECT_TRUE(RunCommand(request).succeeded);
}

TEST(UmapEmbeddingCommandTest, SameSeedAndSettingsProduceIdenticalOutput)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_determinism" };
    const auto input_path{ temp_dir.path() / "local_fitting_result_repeat.csv" };
    WriteCsv(input_path, MakeRows(8));

    auto first_request{ MakeRequest(input_path, temp_dir.path() / "first") };
    auto second_request{ MakeRequest(input_path, temp_dir.path() / "second") };
    ASSERT_TRUE(RunCommand(first_request).succeeded);
    ASSERT_TRUE(RunCommand(second_request).succeeded);

    const auto first_output{ ReadLines(
        first_request.output_dir / "umap_embedding_repeat.csv") };
    const auto second_output{ ReadLines(
        second_request.output_dir / "umap_embedding_repeat.csv") };
    EXPECT_EQ(first_output, second_output);
}

TEST(UmapEmbeddingCommandTest, ZScoreMakesPositiveAffineFeatureChangesInvariant)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_affine" };
    const auto base_path{ temp_dir.path() / "base.csv" };
    const auto transformed_path{ temp_dir.path() / "transformed.csv" };
    WriteCsv(base_path, MakeRows(8));
    WriteCsv(transformed_path, MakeRows(8, std::nullopt, true));

    auto base_request{ MakeRequest(base_path, temp_dir.path() / "base_output") };
    auto transformed_request{
        MakeRequest(transformed_path, temp_dir.path() / "transformed_output")
    };
    ASSERT_TRUE(RunCommand(base_request).succeeded);
    ASSERT_TRUE(RunCommand(transformed_request).succeeded);

    const auto base_coordinates{
        ReadCoordinates(base_request.output_dir / "umap_embedding_base.csv")
    };
    const auto transformed_coordinates{
        ReadCoordinates(transformed_request.output_dir / "umap_embedding_transformed.csv")
    };
    ASSERT_EQ(base_coordinates.size(), transformed_coordinates.size());
    for (std::size_t index = 0; index < base_coordinates.size(); ++index)
    {
        EXPECT_NEAR(base_coordinates[index][0], transformed_coordinates[index][0], 1e-12);
        EXPECT_NEAR(base_coordinates[index][1], transformed_coordinates[index][1], 1e-12);
    }
}

TEST(UmapEmbeddingCommandTest, ReportsOutputWriteFailuresOnFolderOption)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_write_failure" };
    const auto input_path{ temp_dir.path() / "local_fitting_result_blocked.csv" };
    const auto output_dir{ temp_dir.path() / "output" };
    WriteCsv(input_path, MakeRows(6));
    std::filesystem::create_directories(output_dir / "umap_embedding_blocked.csv");

    const auto result{ RunCommand(MakeRequest(input_path, output_dir)) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasIssue(result, "-o,--folder", "Failed to open UMAP output file"));
}

TEST(UmapEmbeddingCommandTest, CliFlagsUseTheSameExecutionPath)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_cli" };
    const auto input_path{ temp_dir.path() / "local_fitting_result_cli.csv" };
    const auto output_dir{ temp_dir.path() / "output" };
    WriteCsv(input_path, MakeRows(6));

    std::vector<std::string> arguments{
        "RHBM-GEM",
        "umap_embedding",
        "--input", input_path.string(),
        "--folder", output_dir.string(),
        "--neighbors", "3",
        "--min-dist", "0.2",
        "--epochs", "5",
        "--seed", "99",
        "--jobs", "1",
    };
    std::vector<char *> argv;
    argv.reserve(arguments.size());
    for (auto & argument : arguments) argv.push_back(argument.data());

    EXPECT_EQ(RunCommandCLI(static_cast<int>(argv.size()), argv.data()), 0);
    EXPECT_TRUE(std::filesystem::exists(output_dir / "umap_embedding_cli.csv"));
}
