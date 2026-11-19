#include "detail/CommandRunner.hpp"
#include "core/classifier/UmapEmbeddingAlgorithm.hpp"
#include "core/classifier/UmapEmbeddingFeatures.hpp"
#include "core/classifier/UmapEmbeddingOutput.hpp"
#include "core/classifier/UmapEmbeddingPlot.hpp"

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace rhbm_gem::core {

namespace {

using classifier::PreparedUmapInput;
using classifier::UmapEmbeddingResult;

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

void NormalizeAndValidateRequest(
    CommandRunner<UmapEmbeddingRequest> & runner,
    UmapEmbeddingRequest & request)
{
    runner.RequireExistingPath(request, &UmapEmbeddingRequest::database_path);
    runner.RequireNonEmptyList(request, &UmapEmbeddingRequest::model_key_tag);

    std::error_code file_error;
    if (!request.database_path.empty()
        && std::filesystem::exists(request.database_path, file_error)
        && !std::filesystem::is_regular_file(request.database_path, file_error))
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::database_path,
            "database_path must name a regular file: " + request.database_path.string());
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

void ValidatePreparedRequest(
    CommandRunner<UmapEmbeddingRequest> & runner,
    const UmapEmbeddingRequest & request,
    std::optional<PreparedUmapInput> & prepared_input,
    std::filesystem::path & output_path)
{
    std::unique_ptr<ModelObject> model_object;
    try
    {
        DataRepository repository{ request.database_path };
        try
        {
            model_object = repository.LoadModel(request.model_key_tag);
        }
        catch (const std::exception & error)
        {
            runner.AddFieldValidationError(
                &UmapEmbeddingRequest::model_key_tag,
                "Failed to load model key '" + request.model_key_tag
                    + "': " + error.what());
            return;
        }
    }
    catch (const std::exception & error)
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::database_path,
            "Failed to open or validate database '"
                + request.database_path.string() + "': " + error.what());
        return;
    }

    std::vector<classifier::LocalFittingFeatureRow> source_rows;
    try
    {
        source_rows = classifier::BuildLocalFittingFeatureRows(*model_object);
    }
    catch (const std::exception & error)
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::model_key_tag,
            "Failed to build local fitting features for model key '"
                + request.model_key_tag + "': " + error.what());
        return;
    }

    std::string preparation_error;
    prepared_input = classifier::BuildAndStandardizeInput(
        *model_object, std::move(source_rows), request.model_key_tag,
        request.num_neighbors, preparation_error);
    if (!prepared_input)
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::model_key_tag,
            preparation_error);
        return;
    }

    output_path = classifier::GetEmbeddingOutputPath(request.output_dir, request.model_key_tag);

    if (!prepared_input->constant_features.empty())
    {
        runner.AddFieldNormalizationWarning(
            &UmapEmbeddingRequest::model_key_tag,
            "Constant feature columns were standardized to zero: "
                + JoinFeatureNames(prepared_input->constant_features) + ".");
    }
    if (prepared_input->effective_neighbors != request.num_neighbors)
    {
        runner.AddFieldNormalizationWarning(
            &UmapEmbeddingRequest::num_neighbors,
            "num_neighbors was limited to "
                + std::to_string(prepared_input->effective_neighbors)
                + " for " + std::to_string(prepared_input->feature_rows.size())
                + " data rows.");
    }
}

bool ExecutePreparedRequest(
    CommandRunner<UmapEmbeddingRequest> & runner,
    const UmapEmbeddingRequest & request,
    const PreparedUmapInput & prepared,
    const std::filesystem::path & output_path)
{
    UmapEmbeddingResult first;
    try
    {
        first.coordinates = classifier::RunUmap(prepared.features, prepared.effective_neighbors,
            request.min_dist, request.num_epochs, request.random_seed, request.job_count);
    }
    catch (const std::exception & error)
    {
        runner.RequirePrepareCondition(false,
            "UMAP embedding failed for model key '" + request.model_key_tag
                + "' in database '" + request.database_path.string()
                + "': " + error.what());
        return false;
    }

    if (!classifier::HasFiniteCoordinates(first.coordinates))
    {
        runner.RequirePrepareCondition(false,
            "UMAP produced a non-finite embedding coordinate for model key '"
                + request.model_key_tag + "' in database '"
                + request.database_path.string() + "'.");
        return false;
    }

    try
    {
        first.mean_shift_coordinates = classifier::RunUmapMeanShift(
            first.coordinates, prepared.feature_rows, request.job_count);
    }
    catch (const std::exception & error)
    {
        runner.RequirePrepareCondition(false,
            "UMAP mean shift failed for model key '" + request.model_key_tag
                + "' in database '" + request.database_path.string()
                + "': " + error.what());
        return false;
    }
    if (!classifier::HasFiniteCoordinates(first.mean_shift_coordinates))
    {
        runner.RequirePrepareCondition(false,
            "UMAP mean shift produced a non-finite coordinate for model key '"
                + request.model_key_tag + "' in database '"
                + request.database_path.string() + "'.");
        return false;
    }

    UmapEmbeddingResult second;
    try
    {
        first.binning = classifier::BinMeanShiftEmbedding(prepared.feature_rows, first.mean_shift_coordinates);
        const auto second_features{ classifier::BuildSecondUmapFeatures(
            prepared.feature_rows, first.binning.point_bin_numbers) };
        if (second_features.values.empty())
        {
            Logger::Log(LogLevel::Warning,
                "Second UMAP skipped because all selected feature columns are constant.");
        }
        else
        {
            second.coordinates = classifier::RunUmap(second_features, prepared.effective_neighbors,
                request.min_dist, request.num_epochs, request.random_seed, request.job_count);
            if (!classifier::HasFiniteCoordinates(second.coordinates))
            {
                throw std::runtime_error("Second UMAP produced a non-finite embedding coordinate.");
            }
        }
    }
    catch (const std::exception & error)
    {
        runner.RequirePrepareCondition(false,
            "UMAP binning or second embedding failed for model key '" + request.model_key_tag
                + "' in database '" + request.database_path.string() + "': " + error.what());
        return false;
    }

    if (!second.coordinates.empty())
    {
        try
        {
            second.mean_shift_coordinates = classifier::RunUmapMeanShift(
                second.coordinates, prepared.feature_rows, request.job_count);
            if (!classifier::HasFiniteCoordinates(second.mean_shift_coordinates))
            {
                throw std::runtime_error("Mean shift produced a non-finite coordinate.");
            }
            second.binning = classifier::BinMeanShiftEmbedding(
                prepared.feature_rows, second.mean_shift_coordinates);
        }
        catch (const std::exception & error)
        {
            runner.RequirePrepareCondition(false,
                "Second UMAP mean shift or binning failed for model key '" + request.model_key_tag
                    + "' in database '" + request.database_path.string() + "': " + error.what());
            return false;
        }
    }

    std::string write_error;
    if (!classifier::WriteEmbeddingOutputs(output_path, prepared, first, second, write_error)
        || !classifier::WriteEmbeddingPlots(output_path, prepared, first, second, write_error))
    {
        runner.AddFieldValidationError(&CommandRequestBase::output_dir, write_error);
        return false;
    }
    return true;
}

} // namespace

namespace command_internal {

CommandResult ExecuteUmapEmbeddingCommand(const UmapEmbeddingRequest & request)
{
    CommandRunner<UmapEmbeddingRequest> runner;
    std::optional<PreparedUmapInput> prepared_input;
    std::filesystem::path output_path;
    return runner.Run(
        request,
        NormalizeAndValidateRequest,
        [&prepared_input, &output_path](
            auto & validation_runner, const auto & prepared_request)
        {
            ValidatePreparedRequest(validation_runner, prepared_request, prepared_input, output_path);
        },
        [&runner, &prepared_input, &output_path](const auto & prepared_request)
        {
            return ExecutePreparedRequest(runner, prepared_request, *prepared_input, output_path);
        });
}

} // namespace command_internal

} // namespace rhbm_gem::core
