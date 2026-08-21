#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

#include "core/detail/GaussianEstimatorStages.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <fstream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem::core {
namespace {
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };

RHBMExecutionOptions MakeExecutionOptions(const FitOptions & options)
{
    RHBMExecutionOptions execution_options;
    execution_options.quiet_mode = false;
    execution_options.thread_size = options.thread_size;
    return execution_options;
}

rhbm_trainer::RHBMTrainingOptions MakeTrainingOptions(const FitOptions & options)
{
    rhbm_trainer::RHBMTrainingOptions training_options;
    training_options.execution_options = MakeExecutionOptions(options);
    return training_options;
}

float CalculateAdjustedResponse(double sample_response, double distance, const GaussianModel3D & offset_model)
{
    const auto evaluation{ offset_model.EvaluateAtDistance(distance) };
    return static_cast<float>(sample_response - (evaluation.response - evaluation.signal));
}

LocalPotentialSampleList BuildSamplesForZeroOffsetGaussianFit(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model)
{
    LocalPotentialSampleList adjusted_sampling_entries;
    adjusted_sampling_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto response{
            CalculateAdjustedResponse(static_cast<double>(sample.response), distance, model)
        };
        adjusted_sampling_entries.emplace_back(LocalPotentialSample{ response, sample.point });
    }
    return adjusted_sampling_entries;
}

LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    double offset = 0.0)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols).WithOffset(offset)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde).WithOffset(offset)
    };
    return LocalGaussianResult{
        alpha_r,
        GaussianModel3DWithUncertainty{ ols_model, GaussianModel3DUncertainty{} },
        GaussianModel3DWithUncertainty{ mdpde_model, GaussianModel3DUncertainty{} },
        std::nullopt,
        false,
        0.0,
        fit_result
    };
}

std::vector<LocalGaussianResult> DecodeMemberGaussianResults(
    const RHBMGroupEstimationResult & result,
    const std::vector<double> & member_offset_list)
{
    const auto member_count{ static_cast<std::size_t>(result.beta_posterior_matrix.cols()) };
    if (member_offset_list.size() != member_count)
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }
    if (result.capital_sigma_posterior_list.size() != member_count)
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }
    eigen_validation::RequireVectorSize(
        result.outlier_flag_array, result.beta_posterior_matrix.cols(),
        "outlier_flag_array", "Group Gaussian member result count is inconsistent.");
    eigen_validation::RequireVectorSize(
        result.statistical_distance_array, result.beta_posterior_matrix.cols(),
        "statistical_distance_array", "Group Gaussian member result count is inconsistent.");

    std::vector<LocalGaussianResult> member_results;
    member_results.reserve(member_count);
    for (Eigen::Index i = 0; i < result.beta_posterior_matrix.cols(); i++)
    {
        const auto member_index{ static_cast<std::size_t>(i) };
        const auto offset{ member_offset_list.at(member_index) };
        const auto gaussian{
            linearization_service::DecodeParameterVector(
                result.beta_posterior_matrix.col(i),
                result.capital_sigma_posterior_list.at(member_index))
        };
        const auto gaussian_with_offset{
            GaussianModel3DWithUncertainty{
                gaussian.GetModel().WithOffset(offset),
                gaussian.GetStandardDeviationModel()
            }
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian_with_offset,
            gaussian_with_offset,
            gaussian_with_offset,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i)
        });
    }
    return member_results;
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    const std::vector<double> & member_offset_list)
{
    const auto group_offset{ array_helper::ComputeMedian(member_offset_list) };
    const auto prior{
        linearization_service::DecodeParameterVector(result.mu_prior, result.capital_lambda)
    };
    const auto mean{
        linearization_service::DecodeParameterVector(result.mu_mean).WithOffset(group_offset)
    };
    const auto mdpde{
        linearization_service::DecodeParameterVector(result.mu_mdpde).WithOffset(group_offset)
    };
    auto member_results{ DecodeMemberGaussianResults(result, member_offset_list) };
    return GroupGaussianResult{
        alpha_g,
        mean,
        mdpde,
        GaussianModel3DWithUncertainty{
            prior.GetModel().WithOffset(group_offset),
            prior.GetStandardDeviationModel()
        },
        std::move(member_results)
    };
}

void OutputLocalFittingResultTable(
    const ModelObject & model_object,
    bool peeling_applied,
    const std::filesystem::path & output_path)
{
    const auto table{
        model_object.GetAnalysisView().GetLocalFittingResultCsv(peeling_applied)
    };
    std::ofstream output{ output_path, std::ios::out | std::ios::trunc };
    if (!output.is_open())
    {
        throw std::runtime_error(
            "Failed to open local fitting result CSV file: " + output_path.string());
    }
    output << table;
    output.close();
    if (!output)
    {
        throw std::runtime_error(
            "Failed to write local fitting result CSV file: " + output_path.string());
    }
}

} // namespace

void RunFixedOffsetLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto selected_atom_size{ atom_list.size() };
    auto analysis{ model_object.EditAnalysis() };
    analysis.EnsureSelectedAtomLocalPotentials();
    size_t atom_count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");
    }

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto & atom{ *atom_list[i] };
        const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
        LocalPotentialSampleList sample_entries{ local_view.GetSamplingEntries(stage) };
        GaussianModel3D offset_model{ local_view.GetGaussianResult(stage).mdpde.GetModel() };
        auto result{
            EstimateLocalGaussian(
                sample_entries,
                local_view.GetAlphaR(stage),
                options,
                offset_model)
        };
        analysis.ApplyAtomLocalGaussianResult(stage, atom, std::move(result));

        if (!options.quiet_mode)
        {
#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                atom_count++;
                Logger::ProgressPercent(atom_count, selected_atom_size);
            }
        }
    }
}

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");

    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    std::size_t response_count_min{ std::numeric_limits<std::size_t>::max() };
    for (const auto & sample_entries : sample_entries_list)
    {
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                sample_entries, options.distance_min, options.distance_max));
        const auto response_count{ static_cast<std::size_t>(dataset_list.back().y.size()) };
        if (response_count < response_count_min) response_count_min = response_count;
    }
    auto training_options{ MakeTrainingOptions(options) };
    if (!dataset_list.empty())
    {
        if (response_count_min < 2)
        {
            return training_options.alpha_min;
        }
        if (training_options.subset_size > response_count_min)
        {
            training_options.subset_size = response_count_min;
        }
    }
    return rhbm_trainer::CrossValidationAlphaR(dataset_list, training_options).best_alpha;
}

double TrainAlphaG(
    const std::vector<std::vector<LocalGaussianResult>> & member_result_list,
    const FitOptions & options)
{
    std::vector<std::vector<RHBMParameterVector>> beta_group_list;
    beta_group_list.reserve(member_result_list.size());
    for (const auto & member_results : member_result_list)
    {
        std::vector<RHBMParameterVector> beta_list;
        beta_list.reserve(member_results.size());
        for (const auto & member_result : member_results)
        {
            beta_list.emplace_back(
                linearization_service::EncodeGaussianToParameterVector(member_result.mdpde.GetModel()));
        }
        beta_group_list.emplace_back(std::move(beta_list));
    }

    const auto training_options{ MakeTrainingOptions(options) };
    if (beta_group_list.empty())
    {
        return training_options.alpha_min;
    }

    return rhbm_trainer::CrossValidationAlphaG(beta_group_list, training_options).best_alpha;
}

LocalGaussianDesignTemplate BuildLocalGaussianDesignTemplate(
    const LocalPotentialSampleList & sample_entries,
    double range_min,
    double range_max)
{
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "fit range");

    LocalGaussianDesignTemplate design_template;
    design_template.source_sample_count = sample_entries.size();
    design_template.source_sample_index_list.reserve(sample_entries.size());
    design_template.distance_list.reserve(sample_entries.size());
    for (std::size_t sample_index = 0; sample_index < sample_entries.size(); sample_index++)
    {
        const auto distance{
            static_cast<double>(sample_entries.at(sample_index).point.distance)
        };
        if (distance < range_min || distance > range_max) continue;
        design_template.source_sample_index_list.emplace_back(sample_index);
        design_template.distance_list.emplace_back(distance);
    }

    const auto row_count{
        static_cast<Eigen::Index>(design_template.distance_list.size())
    };
    design_template.design_matrix = RHBMDesignMatrix::Zero(row_count, 2);
    for (Eigen::Index row = 0; row < row_count; row++)
    {
        const auto distance{
            design_template.distance_list.at(static_cast<std::size_t>(row))
        };
        design_template.design_matrix(row, 0) = 1.0;
        design_template.design_matrix(row, 1) = -0.5 * distance * distance;
    }
    return design_template;
}

RHBMMemberDataset BuildLocalGaussianPreparedDataset(
    const LocalGaussianDesignTemplate & design_template,
    const std::vector<double> & sample_response_list,
    const GaussianModel3D & offset_model)
{
    numeric_validation::RequireFinite(offset_model.GetOffset(), "offset");
    if (sample_response_list.size() != design_template.source_sample_count ||
        design_template.source_sample_index_list.size() != design_template.distance_list.size() ||
        design_template.design_matrix.rows() != static_cast<Eigen::Index>(design_template.distance_list.size()) ||
        design_template.design_matrix.cols() != 2)
    {
        throw std::invalid_argument("Prepared local Gaussian design inputs are inconsistent.");
    }

    const auto candidate_count{
        static_cast<Eigen::Index>(design_template.distance_list.size())
    };
    RHBMMemberDataset dataset;
    dataset.X = RHBMDesignMatrix::Zero(candidate_count, 2);
    dataset.y = RHBMResponseVector::Zero(candidate_count);
    Eigen::Index retained_count{ 0 };
    for (std::size_t row = 0; row < design_template.distance_list.size(); row++)
    {
        const auto sample_index{ design_template.source_sample_index_list.at(row) };
        if (sample_index >= sample_response_list.size())
        {
            throw std::invalid_argument("Prepared local Gaussian sample index is out of range.");
        }
        const auto adjusted_response{
            static_cast<double>(CalculateAdjustedResponse(
                sample_response_list.at(sample_index),
                design_template.distance_list.at(row),
                offset_model))
        };
        if (adjusted_response <= 0.0) continue;
        numeric_validation::RequireFinite(
            adjusted_response,
            "response",
            "Member dataset contains non-finite value.");
        dataset.X.row(retained_count) = design_template.design_matrix.row(
            static_cast<Eigen::Index>(row));
        dataset.y(retained_count) = std::log(adjusted_response);
        retained_count++;
    }

    if (retained_count == 0)
    {
        dataset.X = RHBMDesignMatrix::Zero(1, 2);
        dataset.y = RHBMResponseVector::Zero(1);
    }
    else
    {
        dataset.X.conservativeResize(retained_count, 2);
        dataset.y.conservativeResize(retained_count);
    }

    return dataset;
}

LocalGaussianResult EstimateLocalGaussianPrepared(
    const LocalGaussianDesignTemplate & design_template,
    const std::vector<double> & sample_response_list,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model)
{
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    auto dataset{
        BuildLocalGaussianPreparedDataset(design_template, sample_response_list, offset_model)
    };
    const auto result{
        rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, MakeExecutionOptions(options))
    };
    return DecodeLocalGaussianResult(alpha_r, result, offset_model.GetOffset());
}

LocalGaussianResult EstimateLocalGaussian(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model)
{
    const auto design_template{
        BuildLocalGaussianDesignTemplate(sample_entries, options.distance_min, options.distance_max)
    };
    std::vector<double> sample_response_list;
    sample_response_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        sample_response_list.emplace_back(static_cast<double>(sample.response));
    }
    return EstimateLocalGaussianPrepared(
        design_template,
        sample_response_list,
        alpha_r,
        options,
        offset_model);
}

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    double alpha_g,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_g, "alpha_g");

    if (sample_entries_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("sample_entries_list and member_result_list sizes are inconsistent.");
    }

    const auto execution_options{ MakeExecutionOptions(options) };
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    std::vector<RHBMBetaEstimateResult> fit_result_list;
    fit_result_list.reserve(member_result_list.size());
    std::vector<double> member_offset_list;
    member_offset_list.reserve(member_result_list.size());
    for (std::size_t i = 0; i < member_result_list.size(); i++)
    {
        const auto & member_result{ member_result_list.at(i) };
        const auto & member_model{ member_result.mdpde.GetModel() };
        const auto sampling_entries{
            BuildSamplesForZeroOffsetGaussianFit(
                sample_entries_list.at(i),
                member_model)
        };
        auto dataset{
            rhbm_helper::BuildMemberDataset(
                sampling_entries,
                options.distance_min,
                options.distance_max)
        };
        fit_result_list.emplace_back(
            rhbm_helper::EstimateBetaMDPDE(
                member_result.alpha_r,
                dataset,
                execution_options));
        dataset_list.emplace_back(std::move(dataset));
        member_offset_list.emplace_back(member_model.GetOffset());
    }
    const auto group_input{ rhbm_helper::BuildGroupInput(dataset_list, fit_result_list) };
    const auto raw_result{ rhbm_helper::EstimateGroup(alpha_g, group_input, execution_options) };
    return DecodeGroupGaussianResult(alpha_g, raw_result, member_offset_list);
}

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys(stage) };

    size_t count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run local alpha training for " + std::to_string(group_key_list.size()) + " groups.");
    }
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(stage, group_key) };
        analysis.EnsureAtomGroupLocalPotentials(stage, group_key);
        std::vector<LocalPotentialSampleList> sample_entries_list;
        sample_entries_list.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            auto sample_entries{ local_view.GetSamplingEntries(stage) };
            if (!local_view.HasEnoughSamplingEntriesInRange(
                    stage,
                    options.distance_min,
                    options.distance_max,
                    kMinimumAlphaRTrainingSampleCount)) continue;
            sample_entries_list.emplace_back(std::move(sample_entries));
        }
        if (!sample_entries_list.empty())
        {
            const auto alpha_r{ TrainAlphaR(sample_entries_list, options) };
            analysis.SetAtomGroupAlphaR(stage, group_key, alpha_r);
        }
        count++;
        if (!options.quiet_mode)
        {
            Logger::ProgressPercent(count, group_key_list.size());
        }
    }
}

namespace {

void RunGroupAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys(stage) };

    std::vector<std::vector<LocalGaussianResult>> member_result_list;
    member_result_list.reserve(group_key_list.size());
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{
            analysis_view.GetAtomObjectList(stage, group_key)
        };
        if (group_atom_list.size() < kMinimumAlphaGTrainingMemberCount) continue;
        if (group_atom_list.front()->IsMainChainAtom() == false) continue;
        analysis.EnsureAtomGroupLocalPotentials(stage, group_key);

        std::vector<LocalGaussianResult> group_member_results;
        group_member_results.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            group_member_results.emplace_back(local_view.GetGaussianResult(stage));
        }
        member_result_list.emplace_back(std::move(group_member_results));
    }

    const auto alpha_g{ TrainAlphaG(member_result_list, options) };
    analysis.InitializeGroupAlpha(stage, alpha_g);
}

} // namespace

void RunGroupPotentialFitting(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    analysis.EnsureSelectedAtomLocalPotentials();
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run atom group fitting.");
    }

    auto group_key_list{ analysis_view.CollectAtomGroupKeys(stage) };
    auto group_key_size{ group_key_list.size() };
    size_t key_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t k = 0; k < group_key_size; k++)
    {
        auto group_key{ group_key_list[k] };
        const auto & atom_list{ analysis_view.GetAtomObjectList(stage, group_key) };
        const auto alpha_g{ analysis_view.GetAtomAlphaG(stage, group_key) };
        std::vector<LocalPotentialSampleList> sample_entries_list;
        std::vector<LocalGaussianResult> member_result_list;
        sample_entries_list.reserve(atom_list.size());
        member_result_list.reserve(atom_list.size());
        for (const auto & atom : atom_list)
        {
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            sample_entries_list.emplace_back(local_view.GetSamplingEntries(stage));
            member_result_list.emplace_back(local_view.GetGaussianResult(stage));
        }
        const auto result{
            EstimateGroupGaussian(sample_entries_list, member_result_list, alpha_g, options)
        };

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            analysis.ApplyAtomGroupGaussianResult(stage, group_key, result);
            key_count++;
            if (!options.quiet_mode)
            {
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options)
{
    model_object.EditAnalysis().InitializeLocalFittingSeedModels();
    
    RunLocalAlphaTraining(model_object, options, FittingStage::First);
    RunFixedOffsetLocalFitting(model_object, options, FittingStage::First);
    RunGroupAlphaTraining(model_object, options, FittingStage::First);
    RunGroupPotentialFitting(model_object, options, FittingStage::First);

    model_object.EditAnalysis().CopyFittingStageState(FittingStage::First, FittingStage::Second);
    const auto peeling_applied{ RunSecondStageLocalFitting(model_object, options) };

    model_object.EditAnalysis().CopyFittingStageState(FittingStage::Second, FittingStage::Third);
    RunLocalAlphaTraining(model_object, options, FittingStage::Third);
    RunFixedOffsetLocalFitting(model_object, options, FittingStage::Third);
    RunGroupAlphaTraining(model_object, options, FittingStage::Third);
    RunGroupPotentialFitting(model_object, options, FittingStage::Third);
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            model_object.GetAnalysisView().GetGroupPriorSpotSummary(FittingStage::Third));
    }
    if (options.local_fitting_result_csv_path.has_value())
    {
        OutputLocalFittingResultTable(model_object, peeling_applied, *options.local_fitting_result_csv_path);
    }
}

} // namespace rhbm_gem::core
