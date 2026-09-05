#include <gtest/gtest.h>

#include <rhbm_gem/core/CommandSystem.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include "core/detail/LocalFittingFeatures.hpp"
#include "io/sqlite/SQLiteWrapper.hpp"
#include "support/CommandTestHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace rhbm_gem::core;

namespace {

namespace rg = rhbm_gem;

constexpr std::size_t kOutputColumnCount{
    rg::core::detail::kLocalFittingColumnCount + 2
};
constexpr double kFirstAmplitudeBase{ 1.123456789012345 };

struct FeatureModelOptions
{
    std::size_t atom_count{ 8 };
    bool constant_features{ false };
    bool exact_affine_values{ false };
    bool affine_transform{ false };
    bool include_analysis{ true };
    bool omit_raw_samples{ false };
    bool omit_peeling_samples{ false };
    bool invalid_signal_ratio{ false };
    bool infinite_amplitude{ false };
    std::string atom_id{ "CA" };
};

double MaybeTransform(double value, bool transform, double scale, double offset)
{
    return transform ? value * scale + offset : value;
}

double FeaturePattern(
    int serial_id,
    std::size_t feature,
    bool constant_features)
{
    const std::size_t observation{
        constant_features ? 0 : static_cast<std::size_t>(serial_id - 1)
    };
    const std::size_t value{
        (observation + 1) * (feature + 2)
            + (observation * observation + 3 * feature) % (feature + 3)
    };
    return static_cast<double>(value);
}

rg::LocalGaussianResult MakeGaussianResult(
    int serial_id,
    std::size_t stage,
    const FeatureModelOptions & options)
{
    double amplitude{ options.exact_affine_values
        ? FeaturePattern(serial_id, stage, options.constant_features)
        : kFirstAmplitudeBase
            + 0.017 * FeaturePattern(serial_id, stage, options.constant_features) };
    double width{ options.exact_affine_values
        ? FeaturePattern(serial_id, stage + 3, options.constant_features)
        : 0.712345678901234
            + 0.003 * FeaturePattern(serial_id, stage + 3, options.constant_features) };
    double offset{ options.exact_affine_values
        ? FeaturePattern(serial_id, stage + 6, options.constant_features)
        : -0.345678901234567
            + 0.007 * FeaturePattern(serial_id, stage + 6, options.constant_features) };
    const double stage_value{ static_cast<double>(stage) };
    const double scale{ options.exact_affine_values ? 2.0 : 2.0 + stage_value };
    const double transform_offset{
        options.exact_affine_values ? 0.0 : 3.0 + stage_value
    };
    amplitude = MaybeTransform(
        amplitude, options.affine_transform, scale, transform_offset);
    width = MaybeTransform(
        width, options.affine_transform, scale, transform_offset);
    offset = MaybeTransform(
        offset, options.affine_transform, scale, transform_offset);
    if (options.infinite_amplitude && serial_id == 1 && stage == 0)
    {
        amplitude = std::numeric_limits<double>::infinity();
    }

    rg::LocalGaussianResult result;
    const rg::GaussianModel3D model{ amplitude, width, offset };
    result.ols = rg::GaussianModel3DWithUncertainty{
        model, rg::GaussianModel3DUncertainty{} };
    result.mdpde = result.ols;
    return result;
}

std::unique_ptr<rg::ModelObject> BuildFeatureModel(
    const FeatureModelOptions & options = {})
{
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    atoms.reserve(options.atom_count);
    for (std::size_t index = 0; index < options.atom_count; ++index)
    {
        const int serial_id{ static_cast<int>(options.atom_count - index) };
        auto atom{ std::make_unique<rg::AtomObject>() };
        atom->SetSerialID(serial_id);
        atom->SetSequenceID(serial_id);
        atom->SetComponentID("ALA");
        atom->SetAtomID(options.atom_id);
        atom->SetChainID("A");
        atom->SetElement(Element::CARBON);
        const double spacing{ options.constant_features ? 3.0 : 0.75 };
        atom->SetPosition(spacing * static_cast<double>(serial_id), 0.0, 0.0);
        atoms.emplace_back(std::move(atom));
    }

    auto model{ std::make_unique<rg::ModelObject>(std::move(atoms)) };
    model->SetKeyTag("source");
    model->SetPdbID("TEST");
    model->SetResolution(2.0);
    model->SetResolutionMethod("test");
    model->SelectAllAtoms();
    if (!options.include_analysis) return model;

    auto analysis{ model->EditAnalysis() };
    analysis.InitializeFromSelection();
    for (const auto & atom : model->GetAtomList())
    {
        const int serial_id{ atom->GetSerialID() };
        const double varying_serial{ options.constant_features
            ? 1.0 : static_cast<double>(serial_id) };
        double signal_ratio{ options.exact_affine_values
            ? FeaturePattern(serial_id, 9, options.constant_features) / 512.0
            : 0.08 + 0.003 * FeaturePattern(
                serial_id, 9, options.constant_features) };
        double tail_ratio{ options.exact_affine_values
            ? FeaturePattern(serial_id, 10, options.constant_features) / 512.0
            : 0.16 + 0.002 * FeaturePattern(
                serial_id, 10, options.constant_features) };
        constexpr double ratio_scale{ 2.0 };
        const double signal_offset{ options.exact_affine_values ? 0.0 : 0.125 };
        const double tail_scale{ options.exact_affine_values ? 2.0 : 1.5 };
        const double tail_offset{ options.exact_affine_values ? 0.0 : 0.0625 };
        signal_ratio = MaybeTransform(
            signal_ratio, options.affine_transform, ratio_scale, signal_offset);
        tail_ratio = MaybeTransform(
            tail_ratio, options.affine_transform, tail_scale, tail_offset);

        const double raw_signal{ options.invalid_signal_ratio && serial_id == 1
            ? 0.0
            : options.exact_affine_values ? 512.0 : 8.0 + varying_serial };
        const double raw_tail{
            options.exact_affine_values ? 512.0 : 11.0 + varying_serial
        };
        if (!(options.omit_raw_samples && serial_id == 1))
        {
            analysis.SetAtomLocalRawSamplingEntries(*atom, {
                LocalPotentialSample{ raw_signal, SamplingPoint{ 0.5 } },
                LocalPotentialSample{ raw_tail, SamplingPoint{ 1.5 } },
            });
        }
        if (!(options.omit_peeling_samples && serial_id == 1))
        {
            analysis.SetAtomLocalPeelingSamplingEntries(*atom, {
                LocalPotentialSample{
                    raw_signal * (1.0 - signal_ratio), SamplingPoint{ 0.5 } },
                LocalPotentialSample{
                    raw_tail * (1.0 - tail_ratio), SamplingPoint{ 1.5 } },
            });
        }
        analysis.SetAtomLocalNeighborCountForPeeling(
            *atom,
            options.constant_features ? 4 : serial_id % 4 + 1);
        analysis.SetAtomLocalGaussianResult(
            rg::FittingStage::First, *atom, MakeGaussianResult(serial_id, 0, options));
        analysis.SetAtomLocalGaussianResult(
            rg::FittingStage::Second, *atom, MakeGaussianResult(serial_id, 1, options));
        analysis.SetAtomLocalGaussianResult(
            rg::FittingStage::Third, *atom, MakeGaussianResult(serial_id, 2, options));
    }
    return model;
}

void SeedFeatureDatabase(
    const std::filesystem::path & database_path,
    const std::string & model_key,
    const FeatureModelOptions & options = {})
{
    rg::DataRepository repository{ database_path };
    repository.SaveModel(*BuildFeatureModel(options), model_key);
}

UmapEmbeddingRequest MakeRequest(
    const std::filesystem::path & database_path,
    const std::filesystem::path & output_dir,
    std::string model_key = "model")
{
    UmapEmbeddingRequest request;
    request.database_path = database_path;
    request.model_key_tag = std::move(model_key);
    request.output_dir = output_dir;
    request.num_neighbors = 4;
    request.num_epochs = 10;
    request.random_seed = 1234;
    request.job_count = 1;
    return request;
}

std::vector<std::string> ReadLines(const std::filesystem::path & path)
{
    std::ifstream input{ path };
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line))
    {
        lines.push_back(std::move(line));
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
            std::stod(fields[rg::core::detail::kLocalFittingColumnCount]),
            std::stod(fields[rg::core::detail::kLocalFittingColumnCount + 1])
        });
    }
    return coordinates;
}

int RunCli(std::vector<std::string> arguments)
{
    std::vector<char *> argv;
    argv.reserve(arguments.size());
    for (auto & argument : arguments)
    {
        argv.push_back(argument.data());
    }
    return RunCommandCLI(static_cast<int>(argv.size()), argv.data());
}

} // namespace

TEST(UmapEmbeddingCommandTest, DefaultsMatchThePublicContract)
{
    const UmapEmbeddingRequest request;

    EXPECT_EQ(request.database_path, GetDefaultDatabasePath());
    EXPECT_TRUE(request.model_key_tag.empty());
    EXPECT_EQ(request.num_neighbors, 15);
    EXPECT_DOUBLE_EQ(request.min_dist, 0.1);
    EXPECT_EQ(request.num_epochs, 0);
    EXPECT_EQ(request.random_seed, 42u);
}

TEST(UmapEmbeddingCommandTest, LoadsSavedAnalysisWithoutLocalFittingCsv)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_sqlite_success" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };
    const auto output_dir{ temp_dir.path() / "output" };
    const std::string model_key{ "model/key weird" };
    SeedFeatureDatabase(database_path, model_key);

    auto request{ MakeRequest(database_path, output_dir, model_key) };
    request.num_neighbors = 15;
    const auto result{ RunCommand(request) };

    ASSERT_TRUE(result.succeeded);
    EXPECT_TRUE(HasIssue(result, "--neighbors", "limited to 7"));
    EXPECT_FALSE(std::filesystem::exists(
        temp_dir.path() / "local_fitting_result_model_key_weird.csv"));

    const auto output_path{ output_dir / "umap_embedding_model_key_weird.csv" };
    const auto output_lines{ ReadLines(output_path) };
    ASSERT_EQ(output_lines.size(), 9u);
    EXPECT_EQ(
        output_lines.front(),
        rg::core::detail::BuildLocalFittingCsvHeader() + ",umap x,umap y");
    for (std::size_t row = 1; row < output_lines.size(); ++row)
    {
        const auto fields{ SplitFields(output_lines[row]) };
        ASSERT_EQ(fields.size(), kOutputColumnCount);
        EXPECT_EQ(std::stoi(fields[0]), static_cast<int>(row));
        EXPECT_EQ(fields[1], "ALA");
        EXPECT_EQ(fields[2], "CA");
        EXPECT_TRUE(std::isfinite(
            std::stod(fields[rg::core::detail::kLocalFittingColumnCount])));
        EXPECT_TRUE(std::isfinite(
            std::stod(fields[rg::core::detail::kLocalFittingColumnCount + 1])));
    }

    const auto first_fields{ SplitFields(output_lines[1]) };
    const double expected_amplitude{ kFirstAmplitudeBase + 0.034 };
    EXPECT_DOUBLE_EQ(std::stod(first_fields[7]), expected_amplitude);
    EXPECT_GT(first_fields[7].size(), 6u);

#ifdef HAVE_ROOT
    EXPECT_TRUE(std::filesystem::is_regular_file(
        output_dir / "umap_embedding_model_key_weird.pdf"));
#else
    EXPECT_FALSE(std::filesystem::exists(
        output_dir / "umap_embedding_model_key_weird.pdf"));
#endif
}

TEST(UmapEmbeddingCommandTest, ClassifiesUnconfiguredSpotsAsOtherInRootPlot)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_other_spots" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };
    FeatureModelOptions options;
    options.atom_id = "CB";
    SeedFeatureDatabase(database_path, "other", options);

    const auto output_dir{ temp_dir.path() / "output" };
    const auto result{ RunCommand(MakeRequest(database_path, output_dir, "other")) };

    ASSERT_TRUE(result.succeeded);
#ifdef HAVE_ROOT
    EXPECT_TRUE(std::filesystem::is_regular_file(
        output_dir / "umap_embedding_other.pdf"));
#else
    EXPECT_FALSE(std::filesystem::exists(
        output_dir / "umap_embedding_other.pdf"));
#endif
}

TEST(UmapEmbeddingCommandTest, AttributesDatabaseAndModelLoadFailures)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_input_failures" };

    const auto missing_database{ temp_dir.path() / "missing.sqlite" };
    auto result{ RunCommand(MakeRequest(
        missing_database, temp_dir.path() / "missing_output")) };
    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasIssue(result, "-d,--database", "does not exist"));
    EXPECT_FALSE(std::filesystem::exists(missing_database));

    result = RunCommand(MakeRequest(
        temp_dir.path(), temp_dir.path() / "directory_output"));
    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasIssue(result, "-d,--database", "regular file"));

    const auto valid_database{ temp_dir.path() / "valid.sqlite" };
    SeedFeatureDatabase(valid_database, "known");
    result = RunCommand(MakeRequest(
        valid_database, temp_dir.path() / "empty_key_output", ""));
    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasIssue(result, "-k,--model-key", "cannot be empty"));

    result = RunCommand(MakeRequest(
        valid_database, temp_dir.path() / "unknown_key_output", "unknown"));
    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasIssue(result, "-k,--model-key", "Failed to load model key"));

    const auto invalid_schema_database{ temp_dir.path() / "invalid_schema.sqlite" };
    {
        rg::DataRepository repository{ invalid_schema_database };
    }
    {
        rg::SQLiteWrapper database{ invalid_schema_database };
        database.Execute("PRAGMA user_version = 999;");
    }
    result = RunCommand(MakeRequest(
        invalid_schema_database,
        temp_dir.path() / "invalid_schema_output",
        "model"));
    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasIssue(result, "-d,--database", "Unsupported SQLite schema"));

    EXPECT_EQ(command_test::CountFilesWithExtension(temp_dir.path(), ".csv"), 0u);
}

TEST(UmapEmbeddingCommandTest, RejectsIncompleteOrNonFiniteSavedAnalysis)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_incomplete" };

    struct InvalidCase
    {
        std::string name;
        FeatureModelOptions options;
        std::string diagnostic;
    };
    std::vector<InvalidCase> cases;

    FeatureModelOptions missing_analysis;
    missing_analysis.include_analysis = false;
    cases.push_back({ "missing_analysis", missing_analysis, "not available" });

    FeatureModelOptions missing_raw;
    missing_raw.omit_raw_samples = true;
    cases.push_back({ "missing_raw", missing_raw, "signal peeling ratio" });

    FeatureModelOptions missing_peeling;
    missing_peeling.omit_peeling_samples = true;
    cases.push_back({ "missing_peeling", missing_peeling, "signal peeling ratio" });

    FeatureModelOptions invalid_ratio;
    invalid_ratio.invalid_signal_ratio = true;
    cases.push_back({ "invalid_ratio", invalid_ratio, "signal peeling ratio" });

    FeatureModelOptions infinite_amplitude;
    infinite_amplitude.infinite_amplitude = true;
    cases.push_back({ "infinite", infinite_amplitude, "amplitude 1st" });

    FeatureModelOptions too_short;
    too_short.atom_count = 2;
    cases.push_back({ "too_short", too_short, "at least 3 data rows" });

    for (const auto & invalid_case : cases)
    {
        SCOPED_TRACE(invalid_case.name);
        const auto case_dir{ temp_dir.path() / invalid_case.name };
        const auto database_path{ case_dir / "analysis.sqlite" };
        SeedFeatureDatabase(database_path, "model", invalid_case.options);

        const auto output_dir{ case_dir / "output" };
        const auto result{
            RunCommand(MakeRequest(database_path, output_dir))
        };

        EXPECT_FALSE(result.succeeded);
        EXPECT_TRUE(HasIssue(
            result, "-k,--model-key", invalid_case.diagnostic));
        EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".csv"), 0u);
    }
}

TEST(UmapEmbeddingCommandTest, RejectsAllConstantSelectedFeatures)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_constant" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };
    FeatureModelOptions options;
    options.constant_features = true;
    SeedFeatureDatabase(database_path, "constant", options);

    const auto output_dir{ temp_dir.path() / "output" };
    const auto result{
        RunCommand(MakeRequest(database_path, output_dir, "constant"))
    };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasIssue(
        result,
        "-k,--model-key",
        "All 12 selected UMAP feature columns are constant"));
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".csv"), 0u);
}

TEST(UmapEmbeddingCommandTest, RejectsInvalidOptions)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_options" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };
    SeedFeatureDatabase(database_path, "model");

    auto request{ MakeRequest(database_path, temp_dir.path() / "neighbors") };
    request.num_neighbors = 1;
    EXPECT_TRUE(HasIssue(RunCommand(request), "--neighbors"));

    request = MakeRequest(database_path, temp_dir.path() / "min_dist_high");
    request.min_dist = 1.01;
    EXPECT_TRUE(HasIssue(RunCommand(request), "--min-dist"));

    request = MakeRequest(database_path, temp_dir.path() / "min_dist_low");
    request.min_dist = -0.01;
    EXPECT_TRUE(HasIssue(RunCommand(request), "--min-dist"));

    request = MakeRequest(database_path, temp_dir.path() / "min_dist_nan");
    request.min_dist = std::numeric_limits<double>::quiet_NaN();
    EXPECT_TRUE(HasIssue(RunCommand(request), "--min-dist"));

    request = MakeRequest(database_path, temp_dir.path() / "epochs");
    request.num_epochs = -1;
    EXPECT_TRUE(HasIssue(RunCommand(request), "--epochs"));
}

TEST(UmapEmbeddingCommandTest, AcceptsDocumentedOptionBoundaries)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_boundaries" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };
    SeedFeatureDatabase(database_path, "model");

    auto request{ MakeRequest(database_path, temp_dir.path() / "minimum") };
    request.num_neighbors = 2;
    request.min_dist = 0;
    request.num_epochs = 1;
    EXPECT_TRUE(RunCommand(request).succeeded);

    request = MakeRequest(database_path, temp_dir.path() / "maximum");
    request.min_dist = 1;
    request.num_epochs = 0;
    EXPECT_TRUE(RunCommand(request).succeeded);
}

TEST(UmapEmbeddingCommandTest, SameSeedAndSettingsProduceIdenticalOutput)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_determinism" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };
    SeedFeatureDatabase(database_path, "repeat");

    const auto first_request{
        MakeRequest(database_path, temp_dir.path() / "first", "repeat")
    };
    const auto second_request{
        MakeRequest(database_path, temp_dir.path() / "second", "repeat")
    };
    ASSERT_TRUE(RunCommand(first_request).succeeded);
    ASSERT_TRUE(RunCommand(second_request).succeeded);

    EXPECT_EQ(
        ReadLines(first_request.output_dir / "umap_embedding_repeat.csv"),
        ReadLines(second_request.output_dir / "umap_embedding_repeat.csv"));
}

TEST(UmapEmbeddingCommandTest, ZScoreMakesPositiveAffineFeatureChangesInvariant)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_affine" };
    const auto base_database{ temp_dir.path() / "base.sqlite" };
    const auto transformed_database{ temp_dir.path() / "transformed.sqlite" };
    FeatureModelOptions base_options;
    base_options.exact_affine_values = true;
    SeedFeatureDatabase(base_database, "model", base_options);
    FeatureModelOptions transformed_options;
    transformed_options.exact_affine_values = true;
    transformed_options.affine_transform = true;
    SeedFeatureDatabase(transformed_database, "model", transformed_options);

    const auto base_request{
        MakeRequest(base_database, temp_dir.path() / "base_output")
    };
    const auto transformed_request{
        MakeRequest(transformed_database, temp_dir.path() / "transformed_output")
    };
    ASSERT_TRUE(RunCommand(base_request).succeeded);
    ASSERT_TRUE(RunCommand(transformed_request).succeeded);

    const auto base_coordinates{
        ReadCoordinates(base_request.output_dir / "umap_embedding_model.csv")
    };
    const auto transformed_coordinates{
        ReadCoordinates(transformed_request.output_dir / "umap_embedding_model.csv")
    };
    ASSERT_EQ(base_coordinates.size(), transformed_coordinates.size());
    for (std::size_t index = 0; index < base_coordinates.size(); ++index)
    {
        EXPECT_NEAR(base_coordinates[index][0], transformed_coordinates[index][0], 1e-9);
        EXPECT_NEAR(base_coordinates[index][1], transformed_coordinates[index][1], 1e-9);
    }
}

TEST(UmapEmbeddingCommandTest, ReportsOutputWriteFailuresOnFolderOption)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_write_failure" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };
    SeedFeatureDatabase(database_path, "blocked");
    const auto output_dir{ temp_dir.path() / "output" };
    std::filesystem::create_directories(
        output_dir / "umap_embedding_blocked.csv");

    const auto result{
        RunCommand(MakeRequest(database_path, output_dir, "blocked"))
    };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasIssue(
        result, "-o,--folder", "Failed to open UMAP output file"));
}

TEST(UmapEmbeddingCommandTest, CliUsesDatabaseAndModelKeyAndRejectsOldInputFlag)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_cli" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };
    const auto output_dir{ temp_dir.path() / "output" };
    SeedFeatureDatabase(database_path, "cli");

    EXPECT_EQ(RunCli({
        "RHBM-GEM",
        "umap_embedding",
        "--database", database_path.string(),
        "--model-key", "cli",
        "--folder", output_dir.string(),
        "--neighbors", "3",
        "--min-dist", "0.2",
        "--epochs", "5",
        "--seed", "99",
        "--jobs", "1",
    }), 0);
    EXPECT_TRUE(std::filesystem::exists(
        output_dir / "umap_embedding_cli.csv"));

    EXPECT_NE(RunCli({
        "RHBM-GEM",
        "umap_embedding",
        "--input", "obsolete.csv",
        "--database", database_path.string(),
        "--model-key", "cli",
        "--folder", (temp_dir.path() / "old_flag_output").string(),
    }), 0);
}
