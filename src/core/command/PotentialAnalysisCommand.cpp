#include "PotentialAnalysisCommand.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>

namespace {
constexpr std::string_view kModelKey{ "model" };
constexpr std::string_view kMapKey{ "map" };
constexpr std::string_view kModelOption{ "--model" };
constexpr std::string_view kMapOption{ "--map" };
constexpr std::string_view kSimResolutionOption{ "--sim-resolution" };
constexpr std::string_view kSaveKeyOption{ "--save-key" };
constexpr std::string_view kSamplingOption{ "--sampling" };
constexpr std::string_view kSamplingMinOption{ "--sampling-min" };
constexpr std::string_view kSamplingMaxOption{ "--sampling-max" };
constexpr std::string_view kSamplingHeightOption{ "--sampling-height" };
constexpr std::string_view kFitMinOption{ "--fit-min" };
constexpr std::string_view kFitMaxOption{ "--fit-max" };
constexpr std::string_view kAlphaROption{ "--alpha-r" };
constexpr std::string_view kAlphaGOption{ "--alpha-g" };
constexpr std::string_view kSamplingRangeIssue{ "--sampling-range" };
constexpr std::string_view kFitRangeIssue{ "--fit-range" };
}

namespace rhbm_gem {

PotentialAnalysisCommand::PotentialAnalysisCommand() :
    CommandWithProfileOptions<
        PotentialAnalysisCommandOptions,
        CommandId::PotentialAnalysis,
        CommonOptionProfile::DatabaseWorkflow>{},
    m_model_key_tag{ kModelKey }, m_map_key_tag{ kMapKey },
    m_map_object{ nullptr }, m_model_object{ nullptr }
{
}

void PotentialAnalysisCommand::ApplyRequest(const PotentialAnalysisRequest & request)
{
    ApplyCommonRequest(request.common);
    SetModelFilePath(request.model_file_path);
    SetMapFilePath(request.map_file_path);
    MutateOptions([&]() { m_options.is_simulation = request.simulation_flag; });
    SetSimulatedMapResolution(request.simulated_map_resolution);
    SetSavedKeyTag(request.saved_key_tag);
    MutateOptions([&]() { m_options.training_report_dir = request.training_report_dir; });
    MutateOptions([&]() { m_options.use_training_alpha = request.training_alpha_flag; });
    MutateOptions([&]() { m_options.is_asymmetry = request.asymmetry_flag; });
    SetSamplingSize(request.sampling_size);
    SetSamplingRangeMinimum(request.sampling_range_min);
    SetSamplingRangeMaximum(request.sampling_range_max);
    SetSamplingHeight(request.sampling_height);
    SetFitRangeMinimum(request.fit_range_min);
    SetFitRangeMaximum(request.fit_range_max);
    SetAlphaR(request.alpha_r);
    SetAlphaG(request.alpha_g);
}

bool PotentialAnalysisCommand::ExecuteImpl()
{
    if (BuildDataObject() == false) return false;
    RunMapObjectPreprocessing();
    RunModelObjectPreprocessing();

    RunAtomMapValueSampling();
    RunAtomGroupClassification();
    if (m_options.use_training_alpha) RunAtomAlphaTraining();
    else RunLocalAtomFitting(m_options.alpha_r);
    RunAtomPotentialFitting();
    RunExperimentalBondWorkflowIfEnabled();
    SaveDataObject();
    return true;
}

void PotentialAnalysisCommand::ValidateOptions()
{
    ResetPrepareIssues(kSimResolutionOption);
    if (m_options.is_simulation && m_options.resolution_simulation <= 0.0)
    {
        AddValidationError(
            kSimResolutionOption,
            "Expected a positive simulated-map resolution when '--simulation true' is selected.");
    }

    ResetPrepareIssues(kSamplingRangeIssue);
    if (m_options.sampling_range_min > m_options.sampling_range_max)
    {
        AddValidationError(kSamplingRangeIssue,
            "Expected --sampling-min <= --sampling-max.");
    }

    ResetPrepareIssues(kFitRangeIssue);
    if (m_options.fit_range_min > m_options.fit_range_max)
    {
        AddValidationError(kFitRangeIssue,
            "Expected --fit-min <= --fit-max.");
    }
}

void PotentialAnalysisCommand::ResetRuntimeState()
{
    m_map_object.reset();
    m_model_object.reset();
}

void PotentialAnalysisCommand::SetSimulatedMapResolution(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.resolution_simulation,
        value,
        kSimResolutionOption,
        0.0,
        "Simulated-map resolution must be a finite non-negative value.");
}

void PotentialAnalysisCommand::SetSavedKeyTag(const std::string & tag)
{
    MutateOptions([&]()
    {
        ResetParseIssues(kSaveKeyOption);
        if (!tag.empty())
        {
            m_options.saved_key_tag = tag;
            return;
        }

        m_options.saved_key_tag = "model";
        AddValidationError(
            kSaveKeyOption,
            "Saved key tag cannot be empty.",
            ValidationPhase::Parse);
    });
}

void PotentialAnalysisCommand::SetSamplingSize(int value)
{
    SetNormalizedScalarOption(
        m_options.sampling_size,
        value,
        kSamplingOption,
        [](int candidate) { return candidate > 0; },
        1500,
        "Sampling size must be positive, reset to default value = 1500");
}

void PotentialAnalysisCommand::SetSamplingRangeMinimum(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.sampling_range_min,
        value,
        kSamplingMinOption,
        0.0,
        "Minimum sampling range must be a finite non-negative value.");
}

void PotentialAnalysisCommand::SetSamplingRangeMaximum(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.sampling_range_max,
        value,
        kSamplingMaxOption,
        1.5,
        "Maximum sampling range must be a finite non-negative value.");
}

void PotentialAnalysisCommand::SetSamplingHeight(double value)
{
    SetFinitePositiveScalarOption(
        m_options.sampling_height,
        value,
        kSamplingHeightOption,
        0.1,
        "Sampling height must be a finite positive value.");
}

void PotentialAnalysisCommand::SetFitRangeMinimum(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.fit_range_min,
        value,
        kFitMinOption,
        0.0,
        "Minimum fitting range must be a finite non-negative value.");
}

void PotentialAnalysisCommand::SetFitRangeMaximum(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.fit_range_max,
        value,
        kFitMaxOption,
        1.0,
        "Maximum fitting range must be a finite non-negative value.");
}

void PotentialAnalysisCommand::SetAlphaR(double value)
{
    SetFinitePositiveScalarOption(
        m_options.alpha_r,
        value,
        kAlphaROption,
        0.1,
        "Alpha-R must be a finite positive value.");
}

void PotentialAnalysisCommand::SetAlphaG(double value)
{
    SetFinitePositiveScalarOption(
        m_options.alpha_g,
        value,
        kAlphaGOption,
        0.2,
        "Alpha-G must be a finite positive value.");
}

void PotentialAnalysisCommand::SetModelFilePath(const std::filesystem::path & path)
{
    SetRequiredExistingPathOption(m_options.model_file_path, path, kModelOption, "Model file");
}

void PotentialAnalysisCommand::SetMapFilePath(const std::filesystem::path & path)
{
    SetRequiredExistingPathOption(m_options.map_file_path, path, kMapOption, "Map file");
}

} // namespace rhbm_gem
