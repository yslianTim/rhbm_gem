#include <gtest/gtest.h>
#include <boost/json.hpp>

#include <rhbm_gem/core/CommandSystem.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include "core/estimator/joint_component/JointComponentEstimator.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include "data/detail/JointStageAdapter.hpp"
#include "core/estimator/joint_component/Numerics.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include "core/estimator/joint_component/HaloMedianReconstruction.hpp"
#include "support/JointPostprocessingTestHelpers.hpp"
#include "core/classifier/LocalFittingFeatures.hpp"
#include "core/classifier/UmapEmbeddingAlgorithm.hpp"
#include "io/sqlite/SQLiteWrapper.hpp"
#include "support/CommandTestHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace rhbm_gem::core;

namespace {

namespace rg = rhbm_gem;

constexpr std::size_t kOutputColumnCount{
    rg::core::classifier::kLocalFittingColumnCount + 7
};
constexpr std::size_t kBinNumberColumn{ rg::core::classifier::kLocalFittingColumnCount + 4 };
constexpr double kAmplitudeBase{ 1.123456789012345 };

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
    bool infinite_first_amplitude{ false };
    bool vary_joint_widths{ false };
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
        : kAmplitudeBase
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
    if (serial_id == 1
        && options.infinite_first_amplitude && stage == 0)
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

double FeatureJointWidth(int serial_id, bool vary_joint_widths)
{
    return vary_joint_widths ? 0.35 + 0.05 * (serial_id % 7) : 0.5;
}

JointInitialization MakeFeatureInitialization(
    const JointProblemInput & input,
    std::size_t donor_count,
    bool vary_joint_widths = false)
{
    JointInitialization initialization;
    initialization.valid = true;
    initialization.reason = "valid-widths";
    initialization.tau.resize(input.atom_ids.size());
    for (std::size_t index = 0; index < input.atom_ids.size(); ++index)
    {
        const double width{ FeatureJointWidth(std::stoi(input.atom_ids[index]), vary_joint_widths) };
        initialization.tau[index] = width;
        JointInitializationAtom atom;
        atom.id = input.atom_ids[index];
        atom.role=rhbm_gem::FittingRole::OutputTarget;
        atom.mdpde[1] = width;
        atom.ols[1] = width;
        atom.alpha = 0.5;
        atom.sample_count = 12;
        const bool donor = index < donor_count;
        atom.reason = donor ? "valid-width" : "median-fallback";
        if (donor) atom.original_tau = width;
        atom.used_tau = width;
        atom.seed_source = donor ? "fitted" : "median-fallback";
        atom.donor_count = donor ? 0 : donor_count;
        atom.initial_stage_fit_qualified = donor;
        atom.initial_stage_qualification = donor ? rg::RHBMSolveQualification::NativeSuccess :
            rg::RHBMSolveQualification::Unqualified;
        atom.native_status = static_cast<int>(donor ? rg::RHBMEstimationStatus::SUCCESS :
            rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED);
        initialization.atoms.push_back(std::move(atom));
    }
    return initialization;
}

std::unique_ptr<rg::ModelObject> BuildFeatureModel(
    const FeatureModelOptions & options = {})
{
    // Keep extra width donors outside the geometry of the feature rows under test.
    const auto count = options.include_analysis ? std::max<std::size_t>(20, options.atom_count) : options.atom_count;
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    atoms.reserve(count);
    for (std::size_t index = 0; index < count; ++index)
    {
        const int serial_id{ static_cast<int>(count - index) };
        auto atom{ std::make_unique<rg::AtomObject>() };
        atom->SetSerialID(serial_id);
        atom->SetSequenceID(serial_id);
        atom->SetComponentID("ALA");
        atom->SetAtomID(options.atom_id);
        atom->SetChainID("A");
        atom->SetElement(Element::CARBON);
        const double spacing{ options.constant_features ? 3.0 : 0.75 };
        const double shift = static_cast<std::size_t>(serial_id) > options.atom_count ? 100.0 : 0.0;
        atom->SetPosition(shift + spacing * static_cast<double>(serial_id), 0.0, 0.0);
        atoms.emplace_back(std::move(atom));
    }

    auto model{ std::make_unique<rg::ModelObject>(std::move(atoms)) };
    model->SetKeyTag("source");
    model->SetPdbID("TEST");
    model->SetResolution(2.0);
    model->SetResolutionMethod("test");
    model->SelectAllAtoms();
    if (!options.include_analysis) return model;

    JointProblemInput input;
    input.selection_domain.emplace();
    input.selection_domain->contributor_policy = std::string(rg::BufferExpansionContributorPolicy);
    for (const auto & atom : model->GetAtomList())
    {
        input.selection_domain->output_indices.push_back(input.atom_ids.size());
        input.atom_ids.push_back(std::to_string(atom->GetSerialID()));
        auto & support = input.support.emplace_back();
        const double amplitude = MakeGaussianResult(atom->GetSerialID(), 1, options).mdpde.GetModel().GetAmplitude();
        const double width{ FeatureJointWidth(atom->GetSerialID(), options.vary_joint_widths) };
        for (int k = 0; k < 12; ++k)
        {
            const double distance = k == 0 ? 0.5 : k == 1 ? 1.5 : 0.125 * (k - 2);
            const auto kernel = joint_component::EvaluateKernel(distance * distance, width, 2.5);
            const auto row = input.row_ids.size();
            input.row_ids.push_back(std::to_string(row));
            input.observations.push_back(amplitude * kernel.gaussian - 0.1 * kernel.charge);
            support.push_back({row, distance * distance});
        }
    }
    const auto row_count = input.row_ids.size();
    input.selection_domain->buffer_expansion =
        rg::JointBufferExpansionProvenance{row_count, row_count, 0, 0, 0, 0};
    const auto initialization = MakeFeatureInitialization(input, count, options.vary_joint_widths);
    const JointProblem problem(std::move(input));
    auto snapshot = CaptureJointAnalysisResult(FitJointComponentsInitialized(problem, initialization));
    if (!IsQualifiedJointAnalysisResult(snapshot))
        throw std::runtime_error("UMAP feature fixture requires a qualified Joint fit.");

    auto analysis{ model->EditAnalysis() };
    analysis.InitializeFromSelection();
    analysis.ApplyJointResult(std::move(snapshot), "umap-feature-test");
    std::size_t row{};
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
            analysis.SetAtomLocalSamplingEntries(*atom, {
                LocalPotentialSample{ raw_signal, SamplingPoint{ 0.5 } },
                LocalPotentialSample{ raw_tail, SamplingPoint{ 1.5 } },
            });
        }
        rg::PostFitPeelingResult peeling;
        peeling.source = rg::AtomLocalPotentialView::For(*atom).GetJointStageEstimate().source;
        using joint_component_test::postprocessing::MakeGridObservation;
        peeling.samples = {
            {raw_signal * (1.0 - signal_ratio), {}, MakeGridObservation(row, row, .5, raw_signal, atom->GetPosition())},
            {raw_tail * (1.0 - tail_ratio), {}, MakeGridObservation(row + 1, row + 1, 1.5, raw_tail, atom->GetPosition())}};
        if (options.omit_peeling_samples && serial_id == 1)
        {
            peeling.samples.back().response.reset();
            peeling.samples.back().reason = "outside-joint-domain";
        }
        analysis.SetAtomPostFitPeeling(*atom, std::move(peeling));
        analysis.SetAtomInitialStageGaussianResult(
            *atom, MakeGaussianResult(serial_id, 0, options));
        row += 12;
    }
    model->SelectAtoms([&](const rg::AtomObject & atom) {
        return static_cast<std::size_t>(atom.GetSerialID()) <= options.atom_count;
    });
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
    const std::filesystem::path & output_path,
    std::size_t coordinate_column = rg::core::classifier::kLocalFittingColumnCount)
{
    const auto lines{ ReadLines(output_path) };
    std::vector<std::array<double, 2>> coordinates;
    for (std::size_t line = 1; line < lines.size(); ++line)
    {
        const auto fields{ SplitFields(lines[line]) };
        if (fields.size() != kOutputColumnCount) return {};
        coordinates.push_back({
            std::stod(fields[coordinate_column]),
            std::stod(fields[coordinate_column + 1])
        });
    }
    return coordinates;
}

std::vector<std::array<double, 2>> ReferenceMeanShift(
    const std::vector<std::array<double, 2>> & original,
    double alpha = 2.0,
    double bandwidth_ratio = 0.5,
    int iterations = 100)
{
    if (original.empty()) return {};
    const auto [min_x, max_x]{ std::minmax_element(original.begin(), original.end(),
        [](const auto & a, const auto & b) { return a[0] < b[0]; }) };
    const auto [min_y, max_y]{ std::minmax_element(original.begin(), original.end(),
        [](const auto & a, const auto & b) { return a[1] < b[1]; }) };
    const double radius{ std::hypot((*max_x)[0] - (*min_x)[0], (*max_y)[1] - (*min_y)[1]) / 10.0 };
    if (radius == 0.0) return original;
    auto queries{ original };
    for (int iteration = 0; iteration < iterations; ++iteration)
    {
        for (auto & query : queries)
        {
            std::array<double, 2> weighted_sum{};
            double weight_sum{};
            // Scanning in row order also fixes the order of the floating-point sums.
            for (const auto & source : original)
            {
                const double distance{ std::hypot(source[0] - query[0], source[1] - query[1]) };
                if (distance > radius) continue;
                const double distance_in_bandwidths{ (distance / radius) / bandwidth_ratio };
                const double weight{ std::exp(-0.5 * alpha * distance_in_bandwidths * distance_in_bandwidths) };
                weight_sum += weight;
                for (std::size_t dimension = 0; dimension < query.size(); ++dimension)
                {
                    weighted_sum[dimension] += source[dimension] * weight;
                }
            }
            if (weight_sum > 0)
            {
                query[0] = weighted_sum[0] / weight_sum;
                query[1] = weighted_sum[1] / weight_sum;
            }
        }
    }
    return queries;
}

std::map<std::array<long long, 2>, std::vector<std::size_t>> ReferenceBins(
    const std::vector<std::array<double, 2>> & original,
    const std::vector<std::array<double, 2>> & shifted)
{
    auto lower{ original.front() };
    auto upper{ lower };
    for (const auto & [x, y] : original)
    {
        lower = { std::min(lower[0], x), std::min(lower[1], y) };
        upper = { std::max(upper[0], x), std::max(upper[1], y) };
    }
    const auto width{ std::hypot(upper[0] - lower[0], upper[1] - lower[1]) / 20.0 };
    std::map<std::array<long long, 2>, std::vector<std::size_t>> bins;
    for (std::size_t point = 0; point < shifted.size(); ++point)
    {
        const auto & [x, y]{ shifted[point] };
        bins[{
            static_cast<long long>(std::floor((x - lower[0]) / width)),
            static_cast<long long>(std::floor((y - lower[1]) / width))
        }].push_back(point);
    }
    return bins;
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
    FeatureModelOptions options;
    options.infinite_first_amplitude = true;
    SeedFeatureDatabase(database_path, model_key, options);

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
        "serial id,residue,spot,"
        "neighbor distance sum in 2A,distance to closest neighbor,"
        "signal peeling ratio,tail peeling ratio,"
        "amplitude joint,width joint,offset joint,umap x,umap y,mean shift x,mean shift y,"
        "bin number,second umap x,second umap y");
    EXPECT_EQ(kOutputColumnCount, 17u);
    const auto metadata_lines{ ReadLines(output_path.string() + ".metadata.json") };
    const auto metadata{ boost::json::parse(metadata_lines.at(0)) };
    EXPECT_EQ(metadata.at("standardization").as_string(),
        "population-z-score; constant columns zero");
    EXPECT_EQ(metadata.at("second_umap_features").as_array(), (boost::json::array{
        "amplitude joint median of bin", "width joint median of bin", "offset joint median of bin",
        "amplitude joint", "width joint", "offset joint" }));
    EXPECT_FALSE(metadata.as_object().contains("spatial_feature_normalization"));
    const auto & column_semantics{ metadata.at("column_semantics").as_object() };
    EXPECT_EQ(column_semantics.size(), 3u);
    EXPECT_EQ(column_semantics.at("amplitude joint").as_string(),
        "Joint Gaussian amplitude A");
    EXPECT_EQ(column_semantics.at("width joint").as_string(),
        "Joint shared Gaussian/charge-kernel width tau");
    EXPECT_EQ(column_semantics.at("offset joint").as_string(),
        "Joint charge-kernel coefficient C");
    EXPECT_FALSE(column_semantics.contains("amplitude 2nd"));
    EXPECT_FALSE(column_semantics.contains("width 2nd"));
    EXPECT_FALSE(column_semantics.contains("offset 2nd"));
    EXPECT_TRUE(metadata.at("map_normalization").is_null());
    for (std::size_t row = 1; row < output_lines.size(); ++row)
    {
        const auto fields{ SplitFields(output_lines[row]) };
        ASSERT_EQ(fields.size(), kOutputColumnCount);
        EXPECT_EQ(std::stoi(fields[0]), static_cast<int>(row));
        EXPECT_EQ(fields[1], "ALA");
        EXPECT_EQ(fields[2], "CA");
        EXPECT_GT(std::stoull(fields[kBinNumberColumn]), 0u);
        for (std::size_t column = rg::core::classifier::kLocalFittingColumnCount;
            column < kOutputColumnCount; ++column)
        {
            EXPECT_TRUE(std::isfinite(std::stod(fields[column])));
        }
    }

    const auto first_fields{ SplitFields(output_lines[1]) };
    EXPECT_DOUBLE_EQ(std::stod(first_fields[3]), 2.25);
    rg::DataRepository repository{ database_path };
    const auto saved_model = repository.LoadModel(model_key);
    const auto expected_amplitude = rg::AtomLocalPotentialView::For(*saved_model->FindAtomPtr(1))
        .GetJointStageModel().GetAmplitude();
    EXPECT_DOUBLE_EQ(std::stod(first_fields[7]), expected_amplitude);
    EXPECT_GT(first_fields[7].size(), 6u);

    FeatureModelOptions geometry_options;
    geometry_options.atom_count = 7;
    geometry_options.include_analysis = false;
    auto geometry_model{ BuildFeatureModel(geometry_options) };
    constexpr std::array<std::array<double, 3>, 7> positions{{
        { 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 1.0 },
        { 1.5, 0.0, 0.0 },
        { 0.0, 0.0, 1.51 },
        { 0.0, 2.0, 0.0 },
        { 2.01, 0.0, 0.0 },
        { 10.0, 0.0, 0.0 },
    }};
    for (const auto & atom : geometry_model->GetAtomList())
    {
        atom->SetPosition(positions[static_cast<std::size_t>(atom->GetSerialID() - 1)]);
    }
    geometry_model->EditAnalysis().InitializeFromSelection();
    geometry_model->SelectAtoms([](const rg::AtomObject & atom)
    {
        return atom.GetSerialID() == 1 || atom.GetSerialID() == 7;
    });
    const auto geometry_rows{
        rg::core::classifier::BuildLocalFittingFeatureRows(*geometry_model)
    };
    ASSERT_EQ(geometry_rows.size(), 2u);
    EXPECT_EQ(geometry_rows[0].serial_id, 1);
    EXPECT_DOUBLE_EQ(
        geometry_rows[0].features[0], 1.0 + 1.5 + 1.51 + 2.0);
    EXPECT_EQ(geometry_rows[1].serial_id, 7);
    EXPECT_DOUBLE_EQ(geometry_rows[1].features[0], 0.0);

    for (const auto & atom : geometry_model->GetAtomList())
    {
        if (atom->GetSerialID() == 2 || atom->GetSerialID() == 4)
        {
            atom->SetElement(Element::HYDROGEN);
        }
    }
    const auto hydrogen_rows{
        rg::core::classifier::BuildLocalFittingFeatureRows(*geometry_model)
    };
    ASSERT_EQ(hydrogen_rows.size(), geometry_rows.size());
    EXPECT_DOUBLE_EQ(hydrogen_rows[0].features[0], 3.5);
    EXPECT_DOUBLE_EQ(hydrogen_rows[1].features[0], 0.0);
    for (std::size_t row = 0; row < geometry_rows.size(); ++row)
    {
        EXPECT_EQ(hydrogen_rows[row].serial_id, geometry_rows[row].serial_id);
        EXPECT_EQ(hydrogen_rows[row].residue, geometry_rows[row].residue);
        EXPECT_EQ(hydrogen_rows[row].spot, geometry_rows[row].spot);
        for (std::size_t feature = 0; feature < geometry_rows[row].features.size(); ++feature)
        {
            if (feature <= 1) continue;
            const auto expected{ geometry_rows[row].features[feature] };
            if (std::isnan(expected))
            {
                EXPECT_TRUE(std::isnan(hydrogen_rows[row].features[feature]));
            }
            else
            {
                EXPECT_DOUBLE_EQ(hydrogen_rows[row].features[feature], expected);
            }
        }
    }

#ifdef HAVE_ROOT
    EXPECT_TRUE(std::filesystem::is_regular_file(
        output_dir / "umap_embedding_model_key_weird.pdf"));
    EXPECT_TRUE(std::filesystem::is_regular_file(
        output_dir / "umap_embedding_model_key_weird_mean_shift.pdf"));
#else
    EXPECT_FALSE(std::filesystem::exists(
        output_dir / "umap_embedding_model_key_weird.pdf"));
    EXPECT_FALSE(std::filesystem::exists(
        output_dir / "umap_embedding_model_key_weird_mean_shift.pdf"));
#endif
}

TEST(UmapEmbeddingCommandTest, ClassifiesUnconfiguredSpotsAsOtherInRootPlot)
{
    command_test::ScopedTempDir temp_dir{ "umap_embedding_other_spots" };
    std::vector<std::vector<Element>> element_sets{
        { Element::HELIUM, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
          Element::PHOSPHORUS, Element::SULFUR, Element::UNK },
        { Element::CARBON, Element::SULFUR, Element::IRON, Element::UNK },
        {},
    };
    for (int atomic_number = 2; atomic_number <= 25; ++atomic_number)
    {
        element_sets.back().push_back(static_cast<Element>(atomic_number));
    }
    element_sets.back().push_back(Element::UNK);
    for (std::size_t scenario = 0; scenario < element_sets.size(); ++scenario)
    {
        SCOPED_TRACE(scenario);
        const auto & elements{ element_sets[scenario] };
        const auto model_key{ "other_" + std::to_string(scenario) };
        const auto database_path{ temp_dir.path() / (model_key + ".sqlite") };
        FeatureModelOptions options;
        options.atom_id = "CB";
        options.atom_count = 5 + elements.size();
        auto model{ BuildFeatureModel(options) };
        constexpr std::array<std::string_view, 5> configured_spots{ "C", "CA", "N", "O", "O" };
        // The second scenario has an N spot but only unselected nitrogen atoms.
        const std::array<Element, 5> configured_elements{
            Element::CARBON, Element::CARBON, scenario == 1 ? Element::SULFUR : Element::NITROGEN,
            Element::OXYGEN, Element::OXYGEN
        };
        for (std::size_t index = 0; index < options.atom_count; ++index)
        {
            auto * atom{ model->FindAtomPtr(static_cast<int>(index + 1)) };
            if (index < configured_spots.size())
            {
                atom->SetAtomID(std::string(configured_spots[index]));
                atom->SetElement(configured_elements[index]);
                if (index == 4) atom->SetComponentID("HOH");
            }
            else
            {
                atom->SetElement(elements[index - configured_spots.size()]);
            }
        }
        if (scenario == 1)
        {
            auto * atom{ model->FindAtomPtr(static_cast<int>(options.atom_count + 1)) };
            ASSERT_NE(atom, nullptr);
            atom->SetElement(Element::NITROGEN);
        }
        {
            rg::DataRepository repository{ database_path };
            repository.SaveModel(*model, model_key);
        }
        const auto output_dir{ temp_dir.path() / model_key };
        testing::internal::CaptureStdout();
#ifndef HAVE_ROOT
        testing::internal::CaptureStderr();
#endif
        const auto result{ RunCommand(MakeRequest(database_path, output_dir, model_key)) };
        const auto log{ testing::internal::GetCapturedStdout() };
#ifndef HAVE_ROOT
        const auto warnings{ testing::internal::GetCapturedStderr() };
        EXPECT_NE(warnings.find("UMAP scatter plots skipped because ROOT support is unavailable."),
            std::string::npos);
#endif
        ASSERT_TRUE(result.succeeded);
        const auto csv_path{ output_dir / ("umap_embedding_" + model_key + ".csv") };
        const auto lines{ ReadLines(csv_path) };
        ASSERT_EQ(lines.size(), options.atom_count + 1);
        EXPECT_EQ(lines[0], rg::core::classifier::BuildLocalFittingCsvHeader()
            + ",umap x,umap y,mean shift x,mean shift y,bin number,second umap x,second umap y");
        std::vector<std::vector<std::string>> point_rows;
        std::map<Element, std::size_t> umap_element_counts;
        for (std::size_t row = 1; row < lines.size(); ++row)
        {
            auto fields{ SplitFields(lines[row]) };
            ASSERT_EQ(fields.size(), kOutputColumnCount);
            const auto serial_id{ std::stoi(fields[0]) };
            EXPECT_EQ(serial_id, static_cast<int>(row));
            const auto * atom{ model->FindAtomPtr(serial_id) };
            ASSERT_NE(atom, nullptr);
            ++umap_element_counts[atom->GetElement()];
            for (std::size_t column = rg::core::classifier::kLocalFittingColumnCount;
                column < kOutputColumnCount; ++column)
            {
                EXPECT_TRUE(std::isfinite(std::stod(fields[column])));
            }
            point_rows.push_back(std::move(fields));
        }
        constexpr std::string_view bin_header{
            "bin number,bin x,bin y,point count,amplitude joint median,width joint median,offset joint median,"
            "spot class count,1st dominate spot,1st dominate spot count,1st dominate spot count ratio,"
            "2nd dominate spot,2nd dominate spot count,2nd dominate spot count ratio,"
            "3rd dominate spot,3rd dominate spot count,3rd dominate spot count ratio"
        };
        constexpr std::array<std::string_view, 5> configured_labels{ "C", "CA", "N", "O", "O (Water)" };
        constexpr std::array ratio_elements{ Element::CARBON, Element::NITROGEN, Element::OXYGEN };
        for (const bool second : { false, true })
        {
            const auto suffix{ second ? "_second_umap_mean_shift_bins.csv" : "_mean_shift_bins.csv" };
            SCOPED_TRACE(suffix);
            const auto coordinates{ second ? ReadCoordinates(csv_path, kBinNumberColumn + 1)
                : ReadCoordinates(csv_path) };
            ASSERT_EQ(coordinates.size(), point_rows.size());
            const auto expected_bins{ ReferenceBins(coordinates, ReferenceMeanShift(coordinates)) };
            const auto bin_lines{ ReadLines(output_dir / ("umap_embedding_" + model_key + suffix)) };
            ASSERT_EQ(bin_lines.size(), expected_bins.size() + 1);
            std::string expected_header{ bin_header };
            if (second)
            {
                expected_header += ",local carbon ratio,local nitrogen ratio,local oxygen ratio"
                    ",global carbon ratio,global nitrogen ratio,global oxygen ratio";
            }
            EXPECT_EQ(bin_lines.front(), expected_header);
            const std::string log_prefix{ second
                ? "\nSecond UMAP mean shift bin count: " : "\nUMAP mean shift bin count: " };
            EXPECT_NE(log.find(log_prefix + std::to_string(expected_bins.size()) + '\n'), std::string::npos);
            std::array<double, ratio_elements.size()> global_ratio_sums{};
            std::size_t bin_line{ 1 };
            for (const auto & [key, points] : expected_bins)
            {
                const auto fields{ SplitFields(bin_lines[bin_line]) };
                ASSERT_EQ(fields.size(), second ? 23u : 17u);
                EXPECT_EQ(std::stoull(fields[0]), bin_line);
                if (!second)
                {
                    for (const auto point : points)
                    {
                        EXPECT_EQ(std::stoull(point_rows[point][kBinNumberColumn]), bin_line);
                    }
                }
                ++bin_line;
                EXPECT_EQ(std::stoll(fields[1]), key[0]);
                EXPECT_EQ(std::stoll(fields[2]), key[1]);
                EXPECT_EQ(std::stoull(fields[3]), points.size());
                for (std::size_t feature = 0; feature < 3; ++feature)
                {
                    std::vector<double> values;
                    for (const auto point : points)
                    {
                        const auto & value{ point_rows[point][7 + feature] };
                        if (!value.empty()) values.push_back(std::stod(value));
                    }
                    if (values.empty())
                    {
                        EXPECT_TRUE(fields[4 + feature].empty());
                        continue;
                    }
                    std::sort(values.begin(), values.end());
                    EXPECT_DOUBLE_EQ(std::stod(fields[4 + feature]),
                        (values[(values.size() - 1) / 2] + values[values.size() / 2]) / 2.0);
                }
                std::map<std::string, std::size_t> counts;
                std::map<Element, std::size_t> element_counts;
                for (const auto point : points)
                {
                    const auto index{ static_cast<std::size_t>(std::stoul(point_rows[point][0]) - 1) };
                    const auto element{ index < configured_elements.size()
                        ? configured_elements[index] : elements[index - configured_elements.size()] };
                    ++element_counts[element];
                    if (index < configured_labels.size())
                    {
                        ++counts[std::string(configured_labels[index])];
                        continue;
                    }
                    ++counts[element == Element::UNK ? "Other Unknown"
                        : "Side-Chain " + ChemicalDataHelper::GetLabel(element)];
                }
                std::vector<std::pair<std::string, std::size_t>> ranking(counts.begin(), counts.end());
                std::stable_sort(ranking.begin(), ranking.end(), [](const auto & lhs, const auto & rhs)
                {
                    return lhs.second > rhs.second;
                });
                EXPECT_EQ(std::stoull(fields[7]), ranking.size());
                for (std::size_t rank = 0; rank < 3; ++rank)
                {
                    const auto column{ 8 + 3 * rank };
                    const auto expected{ rank < ranking.size() ? ranking[rank]
                        : std::pair<std::string, std::size_t>{} };
                    EXPECT_EQ(fields[column], expected.first);
                    EXPECT_EQ(std::stoull(fields[column + 1]), expected.second);
                    EXPECT_DOUBLE_EQ(std::stod(fields[column + 2]),
                        static_cast<double>(expected.second) / static_cast<double>(points.size()));
                }
                if (second)
                {
                    const double point_count{ static_cast<double>(points.size()) };
                    for (std::size_t ratio = 0; ratio < ratio_elements.size(); ++ratio)
                    {
                        const auto element{ ratio_elements[ratio] };
                        const double count{ static_cast<double>(element_counts[element]) };
                        const auto umap_count{ umap_element_counts[element] };
                        const double global_ratio{ std::stod(fields[20 + ratio]) };
                        EXPECT_DOUBLE_EQ(std::stod(fields[17 + ratio]), count / point_count);
                        EXPECT_DOUBLE_EQ(global_ratio,
                            umap_count == 0 ? 0.0 : count / static_cast<double>(umap_count));
                        global_ratio_sums[ratio] += global_ratio;
                    }
                }
            }
            if (second)
            {
                for (std::size_t ratio = 0; ratio < ratio_elements.size(); ++ratio)
                {
                    EXPECT_NEAR(global_ratio_sums[ratio],
                        umap_element_counts[ratio_elements[ratio]] == 0 ? 0.0 : 1.0, 1e-12);
                }
            }
        }
        for (const auto suffix : { "", "_mean_shift", "_second_umap", "_second_umap_mean_shift" })
        {
            const auto plot_path{ output_dir / ("umap_embedding_" + model_key + suffix + ".pdf") };
#ifdef HAVE_ROOT
            ASSERT_TRUE(std::filesystem::is_regular_file(plot_path));
            EXPECT_GT(std::filesystem::file_size(plot_path), 0u);
#else
            EXPECT_FALSE(std::filesystem::exists(plot_path));
#endif
        }
    }
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
    cases.push_back({ "missing_analysis", missing_analysis, "at least 3 data rows" });

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
    const auto feature_rows{
        rg::core::classifier::BuildLocalFittingFeatureRows(*BuildFeatureModel(options))
    };
    ASSERT_EQ(feature_rows.size(), options.atom_count);
    for (const auto & row : feature_rows)
    {
        EXPECT_DOUBLE_EQ(row.features[0], 0.0);
    }

    const auto output_dir{ temp_dir.path() / "output" };
    const auto result{
        RunCommand(MakeRequest(database_path, output_dir, "constant"))
    };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasIssue(
        result,
        "-k,--model-key",
        "All 3 selected UMAP feature columns are constant"));
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

TEST(UmapEmbeddingCommandTest, MeanShiftUsesFixedRadiusAndOriginalPoints)
{
    command_test::ScopedTempDir temp_dir{ "umap_mean_shift" };
    for (const bool vary_joint_widths : { false, true })
    {
        SCOPED_TRACE(vary_joint_widths);
        for (const std::size_t point_count : { 8u, 200u, 208u })
        {
            SCOPED_TRACE(point_count);
            const auto case_dir{ temp_dir.path()
                / (std::to_string(point_count) + (vary_joint_widths ? "_varied" : "_equal")) };
            const auto database_path{ case_dir / "analysis.sqlite" };
            FeatureModelOptions options;
            options.atom_count = point_count;
            options.vary_joint_widths = vary_joint_widths;
            SeedFeatureDatabase(database_path, "model", options);
            auto request{ MakeRequest(database_path, case_dir / "output") };
            request.num_epochs = 75;
            request.job_count = point_count > 200 ? 2 : 1;
            ASSERT_TRUE(RunCommand(request).succeeded);

            const auto csv_path{ request.output_dir / "umap_embedding_model.csv" };
            const auto original{ ReadCoordinates(csv_path) };
            const auto shifted{ ReadCoordinates(csv_path, rg::core::classifier::kLocalFittingColumnCount + 2) };
            ASSERT_EQ(original.size(), point_count);
            ASSERT_EQ(shifted.size(), point_count);
            const auto lines{ ReadLines(csv_path) };
            const auto header{ SplitFields(lines.front()) };
            const auto width_column{ static_cast<std::size_t>(std::distance(
                header.begin(), std::find(header.begin(), header.end(), "width joint"))) };
            ASSERT_LT(width_column, header.size());
            std::vector<double> widths;
            for (std::size_t point = 0; point < point_count; ++point)
            {
                const auto fields{ SplitFields(lines[point + 1]) };
                EXPECT_EQ(std::stoi(fields[0]), static_cast<int>(point + 1));
                const double width{ std::stod(fields[width_column]) };
                ASSERT_TRUE(std::isfinite(width));
                ASSERT_GT(width, 0.0);
                widths.push_back(width);
            }
            if (vary_joint_widths)
            {
                EXPECT_GT(*std::max_element(widths.begin(), widths.end())
                    - *std::min_element(widths.begin(), widths.end()), 0.05);
                const auto equal_width_csv{ temp_dir.path() / (std::to_string(point_count) + "_equal")
                    / "output" / "umap_embedding_model.csv" };
                EXPECT_EQ(original, ReadCoordinates(equal_width_csv));
                EXPECT_EQ(shifted, ReadCoordinates(equal_width_csv,
                    rg::core::classifier::kLocalFittingColumnCount + 2));
            }
            else
            {
                EXPECT_NEAR(*std::max_element(widths.begin(), widths.end()),
                    *std::min_element(widths.begin(), widths.end()), 1e-10);
            }
            const auto expected_shift{ ReferenceMeanShift(original) };
            for (std::size_t point = 0; point < point_count; ++point)
            {
                for (std::size_t dimension = 0; dimension < 2; ++dimension)
                {
                    EXPECT_TRUE(std::isfinite(shifted[point][dimension]));
                    EXPECT_NEAR(shifted[point][dimension], expected_shift[point][dimension], 1e-10);
                }
            }

            for (const bool second : { false, true })
            {
                SCOPED_TRACE(second);
                const auto coordinates{ second ? ReadCoordinates(csv_path, kBinNumberColumn + 1) : original };
                ASSERT_EQ(coordinates.size(), point_count);
                const auto shifted_coordinates{ second ? ReferenceMeanShift(coordinates) : shifted };
                const auto expected_bins{ ReferenceBins(coordinates, shifted_coordinates) };
                const auto summary{ ReadLines(csv_path.parent_path() / (second
                    ? "umap_embedding_model_second_umap_mean_shift_bins.csv"
                    : "umap_embedding_model_mean_shift_bins.csv")) };
                ASSERT_EQ(summary.size(), expected_bins.size() + 1);
                std::size_t number{};
                for (const auto & [index, points] : expected_bins)
                {
                    const auto fields{ SplitFields(summary[++number]) };
                    ASSERT_GE(fields.size(), 4u);
                    EXPECT_EQ(std::stoull(fields[0]), number);
                    EXPECT_EQ(std::stoll(fields[1]), index[0]);
                    EXPECT_EQ(std::stoll(fields[2]), index[1]);
                    EXPECT_EQ(std::stoull(fields[3]), points.size());
                    if (!second)
                    {
                        for (const auto point : points)
                            EXPECT_EQ(std::stoull(SplitFields(lines[point + 1])[kBinNumberColumn]), number);
                    }
                }
            }
        }
    }
}

TEST(UmapEmbeddingCommandTest, RadiusMeanShiftMatchesBruteForceAcrossThreadCounts)
{
    using Coordinate = std::array<double, 2>;
    const std::vector<Coordinate> slow{ { -1.06, -0.5 }, { 1.06, -0.5 }, { 30.0, 1.0 } };
    const auto expected_slow{ ReferenceMeanShift(slow) };
    EXPECT_GT(std::abs(expected_slow[0][0] - ReferenceMeanShift(slow, 1.0)[0][0]), 1e-6);
    EXPECT_GT(std::abs(expected_slow[0][0] - ReferenceMeanShift(slow, 2.0, 1.0)[0][0]), 1e-6);
    EXPECT_GT(std::abs(expected_slow[0][0] - ReferenceMeanShift(slow, 2.0, 0.5, 99)[0][0]), 1e-6);

    // A 30-by-40 bounding box makes the first pair exactly one radius apart.
    const std::vector<Coordinate> boundary{ { -3.0, -4.0 }, { 0.0, 0.0 }, { 27.0, 36.0 } };
    const std::vector<Coordinate> sparse{ { -3.0, -4.0 }, { 0.0, 0.0 }, { 2.0, 1.0 }, { 3.0, -2.0 } };
    std::vector<Coordinate> dense(210, { -1.0, -0.25 });
    dense.insert(dense.end(), 50, { 1.0, 0.25 });
    dense.push_back({ 40.0, 0.0 });
    const auto expected_dense{ ReferenceMeanShift(dense) };
    EXPECT_GT(std::abs(expected_dense[0][0] - dense[0][0]), 1e-3);

    const std::vector<Coordinate> identical(7, { -0.5, 0.25 });
    const std::array<const std::vector<Coordinate> *, 5> fixtures{ &slow, &boundary, &sparse, &dense, &identical };
    for (const auto * fixture : fixtures)
    {
        const auto & original{ *fixture };
        SCOPED_TRACE(original.size());
        const auto expected{ ReferenceMeanShift(original) };
        std::vector<double> embedding;
        for (const auto & position : original)
            embedding.insert(embedding.end(), position.begin(), position.end());

        const std::vector<rg::core::classifier::LocalFittingFeatureRow> rows(original.size());
        const auto serial{ rg::core::classifier::RunUmapMeanShift(embedding, rows, 1).coordinates };
        ASSERT_EQ(serial.size(), embedding.size());
        for (std::size_t point = 0; point < original.size(); ++point)
        {
            EXPECT_NEAR(serial[2 * point], expected[point][0], 1e-10);
            EXPECT_NEAR(serial[2 * point + 1], expected[point][1], 1e-10);
        }
        for (const int threads : { 2, 4 })
            EXPECT_EQ(rg::core::classifier::RunUmapMeanShift(embedding, rows, threads).coordinates, serial);
    }
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

    for (const std::string suffix : {
        ".csv", ".csv.metadata.json", "_mean_shift_bins.csv", "_second_umap_mean_shift_bins.csv" })
    {
        const auto file_name{ "umap_embedding_repeat" + suffix };
        EXPECT_EQ(ReadLines(first_request.output_dir / file_name),
            ReadLines(second_request.output_dir / file_name)) << file_name;
    }
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

    for (const std::string suffix : { "mean_shift_bins.csv", "second_umap_mean_shift_bins.csv" })
    {
        const auto bin_output_dir{ temp_dir.path() / (suffix + "_blocked") };
        const auto blocked_bin_path{ bin_output_dir / ("umap_embedding_blocked_" + suffix) };
        std::filesystem::create_directories(blocked_bin_path);
        const auto bin_result{ RunCommand(MakeRequest(database_path, bin_output_dir, "blocked")) };
        EXPECT_FALSE(bin_result.succeeded);
        EXPECT_TRUE(HasIssue(bin_result, "-o,--folder", blocked_bin_path.string()));
        EXPECT_TRUE(std::filesystem::is_directory(blocked_bin_path));
    }
}

#ifdef HAVE_ROOT
TEST(UmapEmbeddingCommandTest, ReportsDerivedPlotWriteFailuresOnFolderOption)
{
    command_test::ScopedTempDir temp_dir{ "umap_mean_shift_write_failure" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };
    SeedFeatureDatabase(database_path, "blocked");
    for (const std::string suffix : { "mean_shift", "second_umap", "second_umap_mean_shift" })
    {
        SCOPED_TRACE(suffix);
        const auto output_dir{ temp_dir.path() / suffix };
        const auto blocked_path{ output_dir / ("umap_embedding_blocked_" + suffix + ".pdf") };
        std::filesystem::create_directories(blocked_path);

        const auto result{ RunCommand(MakeRequest(database_path, output_dir, "blocked")) };

        EXPECT_FALSE(result.succeeded);
        EXPECT_TRUE(HasIssue(result, "-o,--folder", "UMAP plot output path is not a regular file"));
        EXPECT_TRUE(HasIssue(result, "-o,--folder", blocked_path.string()));
        EXPECT_TRUE(std::filesystem::is_directory(blocked_path));
        EXPECT_TRUE(std::filesystem::is_regular_file(output_dir / "umap_embedding_blocked.pdf"));
    }
}
#endif

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

TEST(UmapEmbeddingCommandTest, ExcludesMissingRequiredFeaturesAndAllowsMissingAncillaryFields)
{
    command_test::ScopedTempDir dir{"umap_feature_availability"};
    for(int scenario=0;scenario<3;++scenario)
    {
        FeatureModelOptions options;
        options.omit_raw_samples=scenario==0;
        options.omit_peeling_samples=scenario==1;
        options.invalid_signal_ratio=scenario==2;
        const auto path=dir.path()/(std::to_string(scenario)+".sqlite");
        SeedFeatureDatabase(path,"model",options);
        const auto request=MakeRequest(path,dir.path()/std::to_string(scenario));
        EXPECT_TRUE(RunCommand(request).succeeded);
        EXPECT_TRUE(std::filesystem::exists(request.output_dir/"umap_embedding_model.csv.metadata.json"));
        EXPECT_EQ(ReadLines(request.output_dir/"umap_embedding_model.csv").size(), scenario == 1 ? 8u : 9u);
    }
}

TEST(UmapEmbeddingCommandTest, SavedJointTargetsEmbedWithCoverageExclusionsAndNoInputFiles)
{
    command_test::ScopedTempDir dir{"joint_umap_saved"};
    FeatureModelOptions model_options; model_options.atom_count=22; model_options.include_analysis=false;
    auto model=BuildFeatureModel(model_options);
    JointProblemInput original;
    original.selection_domain.emplace();
    for(std::size_t a=0;a<model->GetAtomList().size();++a)
    {
        original.atom_ids.push_back(std::to_string(model->GetAtomList()[a]->GetSerialID()));
        original.support.emplace_back();
        if(a!=0) original.selection_domain->output_indices.push_back(a);
        for(int k=0;k<12;++k)
        {
            const int serial=model->GetAtomList()[a]->GetSerialID();
            const int x=3*serial+(k==0 ? 2:k==1 ? 6:0);
            const int y=k<2 ? 10:10+k%5+1;
            const int z=k<2 ? 10:10+k/5+1;
            const double dx=.25*(x-3*serial),dy=.25*(y-10),dz=.25*(z-10);
            const double square=dx*dx+dy*dy+dz*dz;
            const auto kernel=rhbm_gem::core::joint_component::EvaluateKernel(square,0.5,2.5);
            original.support.back().push_back({original.observations.size(),square});
            original.row_ids.push_back(std::to_string(x+100*(y+21*z)));
            original.observations.push_back(2*kernel.gaussian-0.1*kernel.charge);
        }
    }
    const auto original_rows=original.row_ids.size();
    const auto fit_rows=original_rows-original.support.front().size();
    const auto excluded_rows=original_rows-fit_rows;
    original.selection_domain->contributor_policy=std::string(rg::BufferExpansionContributorPolicy);
    original.selection_domain->buffer_expansion=
        rg::JointBufferExpansionProvenance{original_rows,fit_rows,excluded_rows,1,1,0};

    JointProblemInput fitting;
    fitting.selection_domain.emplace();
    fitting.selection_domain->contributor_policy=std::string(rg::BufferExpansionContributorPolicy);
    fitting.selection_domain->buffer_expansion=
        rg::JointBufferExpansionProvenance{original_rows,fit_rows,excluded_rows,0,0,0};
    std::vector<std::size_t> original_to_fitting_row(original_rows,fit_rows);
    std::vector<std::size_t> fitting_to_original_rows;
    for(std::size_t atom=1;atom<original.atom_ids.size();++atom)
    {
        fitting.selection_domain->output_indices.push_back(fitting.atom_ids.size());
        fitting.atom_ids.push_back(original.atom_ids[atom]);
        fitting.support.emplace_back();
        for(const auto & support:original.support[atom])
        {
            if(original_to_fitting_row[support.row]==fit_rows)
            {
                original_to_fitting_row[support.row]=fitting.row_ids.size();
                fitting_to_original_rows.push_back(support.row);
                fitting.row_ids.push_back(original.row_ids[support.row]);
                fitting.observations.push_back(original.observations[support.row]);
            }
            fitting.support.back().push_back({original_to_fitting_row[support.row],support.squared_distance});
        }
    }

    const auto initialization = MakeFeatureInitialization(fitting, 20);

    const auto fitting_to_original_contributors=original.selection_domain->output_indices;
    const auto problem=JointProblem(std::move(fitting)).WithReconstructionDomain(
        std::move(original),fitting_to_original_contributors,std::move(fitting_to_original_rows));
    auto snapshot=CaptureJointAnalysisResult(FitJointComponentsInitialized(problem,initialization));
    snapshot.halo_median_reconstruction=rg::core::detail::BuildHaloMedianReconstruction(
        *model,problem,snapshot,"saved-joint-umap").TakeSnapshot();
    snapshot.metadata.grid_size={100,21,21};
    snapshot.metadata.grid_spacing={.25,.25,.25};
    snapshot.metadata.origin={0.,-2.5,-2.5};
    snapshot.metadata.map_normalization=rg::JointMapNormalization{false,false,1};
    const rg::MapObject grid_geometry(snapshot.metadata.grid_size,
        snapshot.metadata.grid_spacing,snapshot.metadata.origin);
    model->EditAnalysis().InitializeFromSelection();
    model->EditAnalysis().ApplyJointStageEstimates(
        rg::data_internal::BuildJointStageEstimates(snapshot, "saved-joint-umap"));
    model->EditAnalysis().SetJointResult(snapshot);
    for(const auto & atom:model->GetAtomList())
    {
        const auto view=rg::AtomLocalPotentialView::For(*atom);
        if(view.GetJointStageEstimate().source.role!=rg::FittingRole::OutputTarget) continue;
        rg::PostFitPeelingResult peeling; peeling.source=view.GetJointStageEstimate().source;
        const auto contributor=static_cast<std::size_t>(std::find(snapshot.atom_ids.begin(),snapshot.atom_ids.end(),
            std::to_string(atom->GetSerialID()))-snapshot.atom_ids.begin());
        const auto & input=problem.OriginalInput();
        const int serial=atom->GetSerialID();
        const std::array raw{8.0 + serial, 11.0 + serial};
        const std::array ratio{0.08 + 0.003 * FeaturePattern(serial, 9, false),
            0.16 + 0.002 * FeaturePattern(serial, 10, false)};
        for(std::size_t i=0;i<raw.size();++i)
        {
            const auto & support=input.support[contributor][i];
            const auto grid=std::stoull(input.row_ids[support.row]);
            peeling.samples.push_back({raw[i] * (1.0 - ratio[i]), {},
                {support.row,grid,support.squared_distance,input.observations[support.row],grid_geometry.GetGridPosition(grid)}});
        }
        if(atom->GetSerialID()==2) {peeling.samples[1].response.reset(); peeling.samples[1].reason="outside-joint-domain";}
        model->EditAnalysis().SetAtomPostFitPeeling(*atom,peeling);
    }
    const auto path=dir.path()/"joint.sqlite";
    {rg::DataRepository repository(path); repository.SaveModel(*model,"model");}
    model.reset();
    const auto request=MakeRequest(path,dir.path()/"output");
    ASSERT_TRUE(RunCommand(request).succeeded);
    const auto csv=ReadLines(request.output_dir/"umap_embedding_model.csv");
    ASSERT_EQ(csv.size(),21); // twenty valid targets + header; one Halo and one coverage exclusion.
    const auto metadata_lines=ReadLines(request.output_dir/"umap_embedding_model.csv.metadata.json");
    const auto metadata=boost::json::parse(metadata_lines.at(0));
    EXPECT_EQ(metadata.at("estimator").as_string(),"joint-components");
    EXPECT_EQ(metadata.at("peeling_mode").as_string(),"original-grid-v1");
    EXPECT_NE(metadata.at("column_semantics").at("offset joint").as_string().find("charge-kernel coefficient C"),
        std::string::npos);
    const auto & map_normalization{ metadata.at("map_normalization").as_object() };
    EXPECT_EQ(map_normalization.size(), 3u);
    EXPECT_FALSE(map_normalization.at("requested").as_bool());
    EXPECT_FALSE(map_normalization.at("applied").as_bool());
    EXPECT_DOUBLE_EQ(map_normalization.at("divisor").as_double(), 1.0);
    EXPECT_EQ(metadata.at("included_rows").as_int64(),20);
    EXPECT_EQ(metadata.at("features").as_array().size(),3);
    EXPECT_FALSE(metadata.at("exclusions").as_object().empty());
    EXPECT_EQ(metadata.at("exclusions").as_object().at(
        "unavailable: tail peeling ratio (tail-coverage-incomplete:outside-joint-domain)").as_int64(),1);
}
