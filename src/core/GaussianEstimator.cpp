#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

#include "core/detail/FitStateView.hpp"
#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/Objective.hpp"
#include "core/detail/TerminalFailure.hpp"
#include "core/detail/IterationProcess.hpp"
#include "core/detail/LocalGaussianPreparation.hpp"
#include "core/detail/CouplingGraph.hpp"
#include "core/detail/PerformanceCounters.hpp"
#include "core/detail/SecondStageInitialization.hpp"
#include "core/detail/TransformedGaussianModel.hpp"
#include "core/detail/SecondStageContext.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/Convergence.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem::core {
namespace {
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };

rhbm_trainer::RHBMTrainingOptions MakeTrainingOptions(const FitOptions & options)
{
    rhbm_trainer::RHBMTrainingOptions training_options;
    training_options.execution_options = detail::MakeExecutionOptions(options);
    return training_options;
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
            detail::CalculateAdjustedResponse(
                static_cast<double>(sample.response),
                distance,
                model)
        };
        adjusted_sampling_entries.emplace_back(LocalPotentialSample{ response, sample.point });
    }
    return adjusted_sampling_entries;
}

std::vector<LocalGaussianResult> DecodeMemberGaussianResults(
    const RHBMGroupEstimationResult & result,
    const std::vector<double> & member_offset_list)
{
    const auto member_count{ static_cast<std::size_t>(result.beta_posterior_matrix.cols()) };
    if (member_offset_list.size() != member_count ||
        result.capital_sigma_posterior_list.size() != member_count)
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
    return GroupGaussianResult{
        alpha_g,
        mean,
        mdpde,
        GaussianModel3DWithUncertainty{
            prior.GetModel().WithOffset(group_offset),
            prior.GetStandardDeviationModel()
        },
        DecodeMemberGaussianResults(result, member_offset_list)
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

void RunRegularPotentialFittingStage(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage)
{
    RunLocalAlphaTraining(model_object, options, stage);
    RunFixedOffsetLocalFitting(model_object, options, stage);
    RunGroupAlphaTraining(model_object, options, stage);
    RunGroupPotentialFitting(model_object, options, stage);
}

using detail::FitState;

void ApplyFitState(
    ModelObject & model_object,
    const detail::SecondStageContext & context,
    const FitState & iteration_state)
{
    const auto model_snapshot{
        detail::BuildSecondStageModelSnapshot(context, iteration_state)
    };

    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t i = 0; i < context.size(); i++)
    {
        auto adjusted_sampling_entries{
            detail::BuildSecondStageAdjustedSamples(context.at(i), model_snapshot)
        };
        analysis.ApplyAtomLocalSecondStageResult(
            *context.at(i).atom,
            iteration_state.at(i),
            std::move(adjusted_sampling_entries));
    }
}

void AppendOffsetSummary(std::ostringstream & stream, const FitState & state)
{
    std::size_t finite_count{ 0 };
    std::vector<double> absolute_offset_list;
    absolute_offset_list.reserve(state.size());
    for (const auto & result : state)
    {
        const auto offset{ result.mdpde.GetModel().GetOffset() };
        if (!std::isfinite(offset)) continue;
        finite_count++;
        absolute_offset_list.emplace_back(std::abs(offset));
    }
    double median_absolute_offset{ 0.0 };
    double percentile_absolute_offset{ 0.0 };
    double maximum_absolute_offset{ 0.0 };
    if (!absolute_offset_list.empty())
    {
        median_absolute_offset = array_helper::ComputeMedian(absolute_offset_list);
        percentile_absolute_offset = array_helper::ComputePercentile(absolute_offset_list, 0.99);
        maximum_absolute_offset = std::ranges::max(absolute_offset_list);
    }
    stream << std::scientific << std::setprecision(2)
        << "; offsets finite = " << finite_count << " of " << state.size()
        << ", |C| median/p99/max = "
        << median_absolute_offset << "/"
        << percentile_absolute_offset << "/"
        << maximum_absolute_offset;
}

void AppendAuditSummary(std::ostringstream & stream, const detail::AuditedState & audited_state)
{
    const auto & objective{ audited_state.objective };
    stream << "; audit best source = ";
    if (audited_state.source_iteration != 0)
    {
        stream << "accepted iteration " << audited_state.source_iteration;
    }
    else
    {
        stream << "initial";
    }
    stream
        << std::scientific << std::setprecision(2)
        << ", fixed audit objective fit/tail-weighted/offset/total = "
        << objective.fit_range_residual_objective << "/"
        << objective.GetTailValidationPenalty() << "/"
        << objective.offset_plausibility_penalty << "/"
        << objective.GetTotalObjective()
        << ", tail raw/weight = " << objective.tail_validation_loss << "/"
        << detail::kTailValidationWeight;
}

void LogTerminalFallback(
    bool quiet_mode,
    const detail::IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Completed local fitting after "
        << iteration_state.accepted_iteration_count
        << " accepted iterations with last validated states retained";
    detail::AppendTerminalSummary(
        warning_message,
        iteration_state.terminal_failure_state.terminal_summary);
    AppendOffsetSummary(warning_message, iteration_state.previous_state);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogConverged(
    bool quiet_mode,
    const algorithm::ParameterChangeStats & transformed_change_stats,
    const detail::IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Converged after " << iteration_state.accepted_iteration_count
        << " iterations with percentile log-peak-height change = "
        << transformed_change_stats.percentile_list.at(detail::kLogPeakHeightChangeIndex)
        << ", percentile log-width change = "
        << transformed_change_stats.percentile_list.at(detail::kLogWidthChangeIndex)
        << ", and percentile offset-to-peak-ratio change = "
        << transformed_change_stats.percentile_list.at(detail::kOffsetToPeakRatioChangeIndex);
    AppendOffsetSummary(message, iteration_state.previous_state);
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogMaximumIterations(bool quiet_mode, const detail::IterationState & iteration_state)
{
    if (quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    detail::AppendTerminalSummary(
        warning_message,
        iteration_state.terminal_failure_state.terminal_summary);
    const auto * audit_state{
        iteration_state.best_audit_state.has_value() ? &*iteration_state.best_audit_state : nullptr
    };
    if (audit_state != nullptr)
    {
        warning_message << "; applying best validated audit state";
        AppendAuditSummary(warning_message, *audit_state);
    }
    else
    {
        warning_message << "; applying latest validated state";
    }
    AppendOffsetSummary(
        warning_message,
        audit_state != nullptr ? audit_state->state : iteration_state.previous_state);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogSecondStageSummary(
    bool quiet_mode,
    const detail::IterationState & iteration_state,
    std::string_view stop_reason,
    bool final_uses_best_audit)
{
    if (quiet_mode) return;

    const auto & best_audit_state{ iteration_state.best_audit_state };
    const auto final_uses_polish{
        final_uses_best_audit && best_audit_state.has_value() ?
            best_audit_state->uses_polish :
            detail::UsesPolish(iteration_state.previous_polish_provenance)
    };
    Logger::FinishProgressLine();
    std::ostringstream message;
    message << "Second-stage local fitting summary: accepted_iterations="
        << iteration_state.accepted_iteration_count << ", best_iteration=";
    if (!best_audit_state.has_value())
    {
        message << "unavailable";
    }
    else if (best_audit_state->source_iteration != 0)
    {
        message << best_audit_state->source_iteration;
    }
    else
    {
        message << "initial";
    }
    message << ", stop_reason=" << stop_reason << ", best_audit_objective=";
    if (best_audit_state.has_value())
    {
        message << std::scientific << std::setprecision(2)
            << best_audit_state->objective.GetTotalObjective();
    }
    else
    {
        message << "unavailable";
    }
    message << ", final_uses_polish=";
    message << (final_uses_polish ? "yes" : "no");
    message << ", final_state_source="
        << (final_uses_best_audit ? "best-audit" : "latest-validated") << ".";
    Logger::Log(LogLevel::Info, message.str());
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
        response_count_min = std::min(response_count_min, response_count);
    }
    auto training_options{ MakeTrainingOptions(options) };
    if (!dataset_list.empty())
    {
        if (response_count_min < 2)
        {
            return training_options.alpha_min;
        }
        training_options.subset_size = std::min(
            training_options.subset_size,
            response_count_min);
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

LocalGaussianResult EstimateLocalGaussian(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model)
{
    const auto design_template{
        detail::BuildLocalGaussianDesignTemplate(
            sample_entries,
            options.distance_min,
            options.distance_max)
    };
    std::vector<double> sample_response_list;
    sample_response_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        sample_response_list.emplace_back(static_cast<double>(sample.response));
    }
    return detail::EstimateLocalGaussianPrepared(
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

    const auto execution_options{ detail::MakeExecutionOptions(options) };
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

bool RunSecondStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    auto context{ detail::BuildSecondStageContext(model_object, options) };
    detail::StoreSecondStageNeighborCounts(model_object, context);
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    FitState initial_state;
    {
        auto initial_state_build_result{ detail::BuildInitialFitState(context) };
        if (initial_state_build_result.failure != detail::SecondStageInitialStateBuildResult::Failure::None)
        {
            if (!options.quiet_mode)
            {
                const auto unselected_seed_failure{
                    initial_state_build_result.failure == detail::SecondStageInitialStateBuildResult::Failure::UnselectedSeedUnavailable
                };
                Logger::Log(LogLevel::Warning,
                    unselected_seed_failure ?
                        "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
                        "is available for every unselected neighbor atom." :
                        "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
                        "is available for every selected atom.");
                Logger::Log(LogLevel::Info,
                    "Second-stage local fitting summary: accepted_iterations=0, "
                    "best_iteration=unavailable, stop_reason=" +
                    std::string(unselected_seed_failure ? "no-valid-unselected-neighbor-seed" : "no-valid-seed") +
                    ", best_audit_objective=unavailable, final_uses_polish=unavailable, "
                    "final_state_source=unavailable.");
            }
            return false;
        }
        detail::LogSecondStageSeedSelections(
            initial_state_build_result.selection_record_list,
            options.quiet_mode);
        detail::LogUnselectedSecondStageSeedSelections(
            initial_state_build_result.unselected_selection_record_list,
            options.quiet_mode);
        initial_state = std::move(initial_state_build_result.state);
    }
    const auto graph_topology{
        detail::BuildSecondStageGraphTopology(context, initial_state, options.quiet_mode)
    };
    detail::LogGraphTopology(graph_topology, options.quiet_mode);
    auto iteration_state{
        detail::BuildIterationState(context, graph_topology, std::move(initial_state), options)
    };
    detail::PerformanceCounters performance_counters{
        options.quiet_mode,
        context,
        iteration_state.solver_workspace_by_key
    };
    if (iteration_state.best_audit_state.has_value())
    {
        performance_counters.RecordFullStateMaterialization();
    }
    detail::LogObjectiveDomain(iteration_state.objective_domain, options.quiet_mode);
    const auto progress_column_widths{ detail::BuildProgressColumnWidths(context.size()) };
    detail::LogProgressHeader(options.quiet_mode, progress_column_widths);

    std::string_view final_stop_reason;
    bool maximum_iterations_reached{ false };
    for (std::size_t iter = 0; iter < detail::kMaximumIterations; iter++)
    {
        if (iteration_state.active_index_list.empty())
        {
            ApplyFitState(model_object, context, iteration_state.previous_state);
            if (iteration_state.terminal_failure_state.HasFailures())
            {
                LogTerminalFallback(options.quiet_mode, iteration_state);
            }
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Info,
                    "Skip 2nd-stage local atom fitting because no atoms are selected.");
            }
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state,
                "terminal-isolation",
                false);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        auto iteration_result{
            detail::RunIteration(
                context,
                graph_topology,
                options,
                iter + 1,
                iteration_state,
                performance_counters)
        };
        if (iteration_result.objective_domain_changed)
        {
            detail::LogObjectiveDomain(iteration_state.objective_domain, options.quiet_mode, true);
        }
        detail::LogAcceptedBacktrackingDiagnostics(options.quiet_mode, iteration_result.diagnostics);
        if (iteration_result.all_rejected_resolution.has_value())
        {
            detail::LogRejectedClusterDiagnostics(
                options.quiet_mode,
                iteration_result.diagnostics.rejected_cluster_diagnostic_list);
            detail::LogIterationProgress(
                options.quiet_mode,
                progress_column_widths,
                iteration_result.progress);
            detail::LogAllRejectedResolution(
                options.quiet_mode,
                iteration_result.diagnostics.trust_region_update,
                *iteration_result.all_rejected_resolution);
            if (*iteration_result.all_rejected_resolution == detail::AllRejectedResolution::Retry)
            {
                continue;
            }

            final_stop_reason = detail::GetAllRejectedResolutionText(*iteration_result.all_rejected_resolution);
            break;
        }

        detail::LogRejectedClusterDiagnostics(
            options.quiet_mode,
            iteration_result.diagnostics.rejected_cluster_diagnostic_list);
        detail::LogIterationProgress(
            options.quiet_mode,
            progress_column_widths,
            iteration_result.progress);

        if (iteration_result.audit_patience_exhausted)
        {
            final_stop_reason = "audit-patience";
            break;
        }

        if (iteration_result.converged)
        {
            ApplyFitState(model_object, context, iteration_state.previous_state);
            if (iteration_state.terminal_failure_state.HasFailures())
            {
                LogTerminalFallback(options.quiet_mode, iteration_state);
            }
            else
            {
                LogConverged(
                    options.quiet_mode,
                    iteration_result.transformed_change_stats,
                    iteration_state);
            }
            LogSecondStageSummary(
                options.quiet_mode,
                iteration_state,
                "converged",
                false);
            RunGroupPotentialFitting(model_object, options, FittingStage::Second);
            return true;
        }

        if (iter + 1 == detail::kMaximumIterations)
        {
            final_stop_reason = "maximum-iterations";
            maximum_iterations_reached = true;
            break;
        }
    }

    if (!final_stop_reason.empty())
    {
        const auto * audit_state{
            iteration_state.best_audit_state.has_value() ? &*iteration_state.best_audit_state : nullptr
        };
        const auto & final_state{
            audit_state != nullptr ? audit_state->state : iteration_state.previous_state
        };
        ApplyFitState(model_object, context, final_state);
        if (maximum_iterations_reached)
        {
            LogMaximumIterations(options.quiet_mode, iteration_state);
        }
        LogSecondStageSummary(
            options.quiet_mode,
            iteration_state,
            final_stop_reason,
            audit_state != nullptr);
        RunGroupPotentialFitting(model_object, options, FittingStage::Second);
        return true;
    }
    return false;
}

void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options)
{
    model_object.EditAnalysis().InitializeLocalFittingSeedModels();
    
    RunRegularPotentialFittingStage(model_object, options, FittingStage::First);

    model_object.EditAnalysis().CopyFittingStageState(FittingStage::First, FittingStage::Second);
    const auto peeling_applied{ RunSecondStageLocalFitting(model_object, options) };

    model_object.EditAnalysis().CopyFittingStageState(FittingStage::Second, FittingStage::Third);
    RunRegularPotentialFittingStage(model_object, options, FittingStage::Third);
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
