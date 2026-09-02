#include "detail/CommandRunner.hpp"

#include <rhbm_gem/utils/domain/Logger.hpp>

#include <umappp/umappp.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace rhbm_gem::core {

namespace {

constexpr std::size_t kIdentifierColumnCount{ 3 };
constexpr std::size_t kFeatureCount{ 11 };
constexpr std::size_t kInputColumnCount{ kIdentifierColumnCount + kFeatureCount };
constexpr std::size_t kOutputDimensionCount{ 2 };
constexpr std::string_view kInputHeader{
    "serial id,residue,spot,neighbor count,peeling ratio,"
    "amplitude 1st,amplitude 2nd,amplitude 3rd,"
    "width 1st,width 2nd,width 3rd,"
    "offset 1st,offset 2nd,offset 3rd"
};
constexpr std::array<std::string_view, kFeatureCount> kFeatureNames{
    "neighbor count",
    "peeling ratio",
    "amplitude 1st",
    "amplitude 2nd",
    "amplitude 3rd",
    "width 1st",
    "width 2nd",
    "width 3rd",
    "offset 1st",
    "offset 2nd",
    "offset 3rd",
};

struct PreparedUmapInput
{
    std::vector<std::string> original_rows;
    std::vector<double> standardized_features;
    std::vector<std::string_view> constant_features;
    int effective_neighbors{ 0 };
    std::filesystem::path output_path;
};

std::vector<std::string_view> SplitCsvRow(std::string_view row)
{
    std::vector<std::string_view> fields;
    std::size_t field_start{ 0 };
    while (true)
    {
        const auto delimiter{ row.find(',', field_start) };
        if (delimiter == std::string_view::npos)
        {
            fields.emplace_back(row.substr(field_start));
            return fields;
        }
        fields.emplace_back(row.substr(field_start, delimiter - field_start));
        field_start = delimiter + 1;
    }
}

bool ParseFiniteDouble(std::string_view text, double & value)
{
    if (text.empty()) return false;

    const char * const begin{ text.data() };
    const char * const end{ text.data() + text.size() };
    const auto parse_result{ std::from_chars(begin, end, value, std::chars_format::general) };
    return parse_result.ec == std::errc{}
        && parse_result.ptr == end
        && std::isfinite(value);
}

std::string FormatRowLocation(
    std::size_t line_number,
    std::string_view column_name)
{
    return "line " + std::to_string(line_number)
        + ", column '" + std::string(column_name) + "'";
}

std::filesystem::path BuildOutputPath(const UmapEmbeddingRequest & request)
{
    constexpr std::string_view input_prefix{ "local_fitting_result_" };
    std::string suffix{ request.input_csv_path.stem().string() };
    if (suffix.starts_with(input_prefix))
    {
        suffix.erase(0, input_prefix.size());
    }
    if (suffix.empty())
    {
        suffix = "result";
    }
    return request.output_dir / ("umap_embedding_" + suffix + ".csv");
}

std::string JoinFeatureNames(const std::vector<std::string_view> & names)
{
    std::string result;
    for (std::size_t index = 0; index < names.size(); ++index)
    {
        if (index != 0) result += ", ";
        result += names[index];
    }
    return result;
}

std::optional<PreparedUmapInput> ReadAndStandardizeInput(
    const UmapEmbeddingRequest & request,
    std::string & error_message)
{
    std::ifstream input{ request.input_csv_path };
    if (!input.is_open())
    {
        error_message = "Failed to open local fitting result CSV '"
            + request.input_csv_path.string() + "'.";
        return std::nullopt;
    }

    std::string header;
    if (!std::getline(input, header))
    {
        error_message = "line 1, column 'header': expected the local fitting result header, "
            "but the file is empty.";
        return std::nullopt;
    }
    if (!header.empty() && header.back() == '\r') header.pop_back();
    if (header != kInputHeader)
    {
        error_message = "line 1, column 'header': local fitting result header does not match "
            "the required 14-column order.";
        return std::nullopt;
    }

    std::vector<std::array<double, kFeatureCount>> raw_features;
    std::vector<std::string> original_rows;
    std::array<long double, kFeatureCount> means{};
    std::array<long double, kFeatureCount> sum_squared_differences{};

    std::string row;
    std::size_t line_number{ 1 };
    while (std::getline(input, row))
    {
        ++line_number;
        if (!row.empty() && row.back() == '\r') row.pop_back();
        if (row.find('"') != std::string::npos)
        {
            error_message = "line " + std::to_string(line_number)
                + ", column 'layout': quoted CSV fields are not supported.";
            return std::nullopt;
        }

        const auto fields{ SplitCsvRow(row) };
        if (fields.size() != kInputColumnCount)
        {
            error_message = "line " + std::to_string(line_number)
                + ", column 'layout': expected 14 columns, received "
                + std::to_string(fields.size()) + ".";
            return std::nullopt;
        }

        std::array<double, kFeatureCount> feature_row{};
        for (std::size_t feature = 0; feature < kFeatureCount; ++feature)
        {
            const auto text{ fields[kIdentifierColumnCount + feature] };
            if (!ParseFiniteDouble(text, feature_row[feature]))
            {
                error_message = FormatRowLocation(line_number, kFeatureNames[feature])
                    + ": expected a finite numeric value, received '"
                    + std::string(text) + "'.";
                return std::nullopt;
            }
        }

        original_rows.push_back(row);
        raw_features.push_back(feature_row);
        const long double count{ static_cast<long double>(raw_features.size()) };
        for (std::size_t feature = 0; feature < kFeatureCount; ++feature)
        {
            const long double value{ static_cast<long double>(feature_row[feature]) };
            const long double delta{ value - means[feature] };
            means[feature] += delta / count;
            const long double updated_delta{ value - means[feature] };
            sum_squared_differences[feature] += delta * updated_delta;
        }
    }
    if (input.bad())
    {
        error_message = "Failed while reading local fitting result CSV '"
            + request.input_csv_path.string() + "'.";
        return std::nullopt;
    }
    if (raw_features.size() < 3)
    {
        error_message = "UMAP embedding requires at least 3 data rows; found "
            + std::to_string(raw_features.size()) + ".";
        return std::nullopt;
    }
    if (raw_features.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        error_message = "UMAP embedding input exceeds the supported number of rows.";
        return std::nullopt;
    }

    std::array<long double, kFeatureCount> standard_deviations{};
    std::vector<std::string_view> constant_features;
    for (std::size_t feature = 0; feature < kFeatureCount; ++feature)
    {
        const long double variance{
            sum_squared_differences[feature]
                / static_cast<long double>(raw_features.size())
        };
        if (!std::isfinite(means[feature]) || !std::isfinite(variance) || variance < 0)
        {
            error_message = "column '" + std::string(kFeatureNames[feature])
                + "' could not be standardized to finite values.";
            return std::nullopt;
        }
        standard_deviations[feature] = std::sqrt(variance);
        if (standard_deviations[feature] == 0)
        {
            constant_features.push_back(kFeatureNames[feature]);
        }
    }
    if (constant_features.size() == kFeatureCount)
    {
        error_message = "All 11 UMAP feature columns are constant; no embedding can be created.";
        return std::nullopt;
    }

    std::vector<double> standardized_features(raw_features.size() * kFeatureCount);
    for (std::size_t observation = 0; observation < raw_features.size(); ++observation)
    {
        for (std::size_t feature = 0; feature < kFeatureCount; ++feature)
        {
            double standardized_value{ 0 };
            if (standard_deviations[feature] != 0)
            {
                const long double value{
                    (static_cast<long double>(raw_features[observation][feature]) - means[feature])
                        / standard_deviations[feature]
                };
                standardized_value = static_cast<double>(value);
            }
            if (!std::isfinite(standardized_value))
            {
                error_message = FormatRowLocation(observation + 2, kFeatureNames[feature])
                    + ": Z-score is not finite.";
                return std::nullopt;
            }
            standardized_features[observation * kFeatureCount + feature] = standardized_value;
        }
    }

    PreparedUmapInput prepared;
    prepared.original_rows = std::move(original_rows);
    prepared.standardized_features = std::move(standardized_features);
    prepared.constant_features = std::move(constant_features);
    prepared.effective_neighbors = std::min(
        request.num_neighbors,
        static_cast<int>(prepared.original_rows.size() - 1));
    prepared.output_path = BuildOutputPath(request);
    return prepared;
}

void NormalizeAndValidateRequest(
    CommandRunner<UmapEmbeddingRequest> & runner,
    UmapEmbeddingRequest & request)
{
    runner.RequireExistingPath(request, &UmapEmbeddingRequest::input_csv_path);

    std::error_code file_error;
    if (!request.input_csv_path.empty()
        && std::filesystem::exists(request.input_csv_path, file_error)
        && !std::filesystem::is_regular_file(request.input_csv_path, file_error))
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::input_csv_path,
            "input_csv_path must name a regular file: " + request.input_csv_path.string());
    }
    if (request.num_neighbors < 2)
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::num_neighbors,
            "num_neighbors must be at least 2.");
    }
    if (!std::isfinite(request.min_dist)
        || request.min_dist < 0
        || request.min_dist > 1)
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::min_dist,
            "min_dist must be a finite value within [0, 1].");
    }
    if (request.num_epochs < 0)
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::num_epochs,
            "num_epochs must be non-negative.");
    }
}

bool WriteEmbedding(
    const PreparedUmapInput & prepared,
    const std::vector<double> & embedding,
    std::string & error_message)
{
    std::ofstream output{ prepared.output_path, std::ios::out | std::ios::trunc };
    if (!output.is_open())
    {
        error_message = "Failed to open UMAP output file '"
            + prepared.output_path.string() + "'.";
        return false;
    }

    output << kInputHeader << ",umap x,umap y\n";
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t observation = 0; observation < prepared.original_rows.size(); ++observation)
    {
        output << prepared.original_rows[observation]
            << ',' << embedding[observation * kOutputDimensionCount]
            << ',' << embedding[observation * kOutputDimensionCount + 1]
            << '\n';
    }
    output.close();
    if (!output)
    {
        error_message = "Failed while writing UMAP output file '"
            + prepared.output_path.string() + "'.";
        return false;
    }
    return true;
}

bool ExecutePreparedRequest(
    CommandRunner<UmapEmbeddingRequest> & runner,
    const UmapEmbeddingRequest & request,
    const PreparedUmapInput & prepared)
{
    std::vector<double> embedding(
        prepared.original_rows.size() * kOutputDimensionCount);
    try
    {
        knncolle::VptreeBuilder<int, double, double> neighbor_builder(
            std::make_shared<knncolle::EuclideanDistance<double, double>>());
        umappp::Options options;
        options.num_neighbors = prepared.effective_neighbors;
        options.min_dist = request.min_dist;
        if (request.num_epochs > 0) options.num_epochs = request.num_epochs;
        const auto seed{
            static_cast<umappp::RngEngine::result_type>(request.random_seed)
        };
        options.initialize_seed = seed;
        options.initialize_spectral_irlba_options.seed = seed;
        options.optimize_seed = seed;
        options.num_threads = request.job_count;
        options.num_threads_spectral = request.job_count;
        options.num_threads_optimize = 1;

        auto status{ umappp::initialize(
            static_cast<int>(kFeatureCount),
            static_cast<int>(prepared.original_rows.size()),
            prepared.standardized_features.data(),
            neighbor_builder,
            kOutputDimensionCount,
            embedding.data(),
            options) };
        status.run(embedding.data());
    }
    catch (const std::exception & error)
    {
        runner.RequirePrepareCondition(false,
            "UMAP embedding failed for '" + request.input_csv_path.string()
                + "': " + error.what());
        return false;
    }

    for (const double coordinate : embedding)
    {
        if (std::isfinite(coordinate)) continue;
        runner.RequirePrepareCondition(false,
            "UMAP produced a non-finite embedding coordinate for '"
                + request.input_csv_path.string() + "'.");
        return false;
    }

    std::string write_error;
    if (!WriteEmbedding(prepared, embedding, write_error))
    {
        runner.AddFieldValidationError(&CommandRequestBase::output_dir, write_error);
        return false;
    }

    Logger::Log(LogLevel::Info,
        "UMAP embedding written to: " + prepared.output_path.string());
    return true;
}

} // namespace

namespace command_internal {

CommandResult ExecuteUmapEmbeddingCommand(const UmapEmbeddingRequest & request)
{
    CommandRunner<UmapEmbeddingRequest> runner;
    std::optional<PreparedUmapInput> prepared_input;
    return runner.Run(
        request,
        NormalizeAndValidateRequest,
        [&prepared_input](
            CommandRunner<UmapEmbeddingRequest> & validation_runner,
            const UmapEmbeddingRequest & prepared_request)
        {
            std::string preparation_error;
            prepared_input = ReadAndStandardizeInput(prepared_request, preparation_error);
            if (!prepared_input)
            {
                validation_runner.AddFieldValidationError(
                    &UmapEmbeddingRequest::input_csv_path,
                    preparation_error);
                return;
            }

            if (!prepared_input->constant_features.empty())
            {
                validation_runner.AddFieldNormalizationWarning(
                    &UmapEmbeddingRequest::input_csv_path,
                    "Constant feature columns were standardized to zero: "
                        + JoinFeatureNames(prepared_input->constant_features) + ".");
            }
            if (prepared_input->effective_neighbors != prepared_request.num_neighbors)
            {
                validation_runner.AddFieldNormalizationWarning(
                    &UmapEmbeddingRequest::num_neighbors,
                    "num_neighbors was limited to "
                        + std::to_string(prepared_input->effective_neighbors)
                        + " for " + std::to_string(prepared_input->original_rows.size())
                        + " data rows.");
            }
        },
        [&runner, &prepared_input](const UmapEmbeddingRequest & prepared_request)
        {
            if (!prepared_input) return false;
            return ExecutePreparedRequest(runner, prepared_request, *prepared_input);
        });
}

} // namespace command_internal

} // namespace rhbm_gem::core
