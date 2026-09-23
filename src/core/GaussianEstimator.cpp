#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/StageSummary.hpp"
#include "core/detail/PostFitPeeling.hpp"
#include "core/detail/JointUncertainty.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "core/detail/gaussian_fit/FittingRanges.hpp"
#include "core/detail/gaussian_fit/GaussianModelOperations.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include "core/detail/gaussian_fit/PreparedLocalGaussianFit.hpp"

#include <algorithm>
#include <set>
#include <chrono>
#include <array>
#include <cstddef>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::core {
namespace {
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };
constexpr std::array<Spot, 5> kLocalMDPDESummarySpotList{
    Spot::C, Spot::CA, Spot::CB, Spot::N, Spot::O
};

struct GaussianModelParameterSamples
{
    std::vector<double> amplitude_list, width_list, offset_list;
    std::size_t unavailable{}, not_converged{};
};
} // namespace

StageProvenance CollectStageProvenance(const std::vector<const AtomObject *> & atoms)
{
    std::set<std::string> methods, modes;
    for (const auto * atom : atoms)
    {
        const auto view = AtomLocalPotentialView::For(*atom);
        if (!view.IsAvailable()) { methods.insert("unknown"); modes.insert("unknown"); continue; }
        const auto method = view.GetStageEstimate(FittingStage::Second).source.method;
        methods.insert(method == EstimateMethod::JointComponents ? "joint-components" :
            method == EstimateMethod::Peeling ? "two-stage" : "unknown");
        if (view.GetPostFitPeeling()) modes.insert(view.GetPostFitPeeling()->mode);
        else modes.insert(method == EstimateMethod::Peeling && !view.GetPeelingSamplingEntries().empty() ? "iterative" : "unknown");
    }
    const auto label = [](const auto & values) -> std::string {
        return values.empty() ? "unknown" : values.size() == 1 ? *values.begin() : "mixed";
    };
    return {label(methods), label(modes)};
}

std::string BuildSecondStageSpotSummary(const ModelObject & model_object)
{
    std::map<Spot, GaussianModelParameterSamples> spots;
    std::vector<const AtomObject *> population;
    for (const auto * atom : model_object.GetSelectedAtoms())
    {
        const auto spot = atom->GetSpot();
        if (std::find(kLocalMDPDESummarySpotList.begin(), kLocalMDPDESummarySpotList.end(), spot) ==
            kLocalMDPDESummarySpotList.end()) continue;
        const auto view = AtomLocalPotentialView::For(*atom);
        const bool joint = view.IsAvailable() && view.GetStageEstimate(FittingStage::Second).source.method == EstimateMethod::JointComponents;
        if (joint && view.GetStageEstimate(FittingStage::Second).source.role != FittingRole::Target) continue;
        population.push_back(atom);
        auto & samples = spots[spot];
        if (!view.HasFinalModel(FittingStage::Second)) { ++samples.unavailable; continue; }
        const auto & estimate = view.GetStageEstimate(FittingStage::Second);
        if (joint && estimate.convergence != JointCheckStatus::Passed) ++samples.not_converged;
        const auto & point = *estimate.point;
        samples.amplitude_list.push_back(point.GetAmplitude());
        samples.width_list.push_back(point.GetWidth());
        samples.offset_list.push_back(point.GetOffset());
    }
    std::ostringstream summary;
    summary << "Second-stage estimate summary by Spot:\nEstimator: "
        << CollectStageProvenance(population).estimator
        << "\nPopulation: selected targets; s.d.: between-atom dispersion"
        << "\n| Spot | valid | not-converged | unavailable | A mean / s.d. | B mean / s.d. | C charge coefficient mean / s.d. |";
    for (const auto & [spot, samples] : spots)
    {
        summary << "\n| " << ChemicalDataHelper::GetLabel(spot) << " | " << samples.amplitude_list.size()
            << " | " << samples.not_converged << " | " << samples.unavailable;
        for (const auto * values : {&samples.amplitude_list, &samples.width_list, &samples.offset_list})
        {
            if (values->empty()) { summary << " | unavailable"; continue; }
            const auto mean = array_helper::ComputeMean(values->data(), values->size());
            summary << " | " << std::fixed << std::setprecision(2) << mean << " / "
                << array_helper::ComputeStandardDeviation(values->data(), values->size(), mean);
        }
        summary << " |";
    }
    if (spots.empty()) summary << "\nNo matching selected targets available.";
    return summary.str();
}

namespace {
rhbm_trainer::RHBMTrainingOptions MakeTrainingOptions(const FitOptions & options)
{
    rhbm_trainer::RHBMTrainingOptions training_options;
    training_options.execution_options = RHBMExecutionOptions{
        .thread_size = options.thread_size
    };
    return training_options;
}

std::vector<double> CollectSampleResponses(const LocalPotentialSampleList & sample_entries)
{
    std::vector<double> response_list;
    response_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        response_list.emplace_back(sample.response);
    }
    return response_list;
}

std::vector<GroupGaussianMemberResult> DecodeMemberGaussianResults(
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

    std::vector<GroupGaussianMemberResult> member_results;
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
            detail::WithPreservedUncertaintyOffset(gaussian, offset)
        };
        member_results.emplace_back(GroupGaussianMemberResult{
            gaussian_with_offset,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i),
            {},
            true,
            std::nullopt
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
        detail::WithPreservedUncertaintyOffset(prior, group_offset),
        DecodeMemberGaussianResults(result, member_offset_list)
    };
}

void RunGroupAlphaTraining(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys() };

    std::vector<std::vector<GaussianModel3D>> model_group_list;
    model_group_list.reserve(group_key_list.size());
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key) };
        if (group_atom_list.size() < kMinimumAlphaGTrainingMemberCount) continue;
        if (group_atom_list.front()->IsMainChainAtom() == false) continue;
        analysis.EnsureAtomGroupLocalPotentials(group_key);

        std::vector<GaussianModel3D> group_member_models;
        group_member_models.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            const auto local_view{ AtomLocalPotentialView::For(*atom) };
            group_member_models.emplace_back(
                local_view.GetEstimateMDPDE(FittingStage::Second));
        }
        model_group_list.emplace_back(std::move(group_member_models));
    }

    const auto alpha_g{ TrainAlphaG(model_group_list, options) };
    analysis.InitializeGroupAlpha(alpha_g);
}

} // namespace

void RunFixedOffsetLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage)
{
    model_object.EditAnalysis().EnsureSelectedAtomLocalPotentials();
    RunFixedOffsetLocalFitting(model_object, options, stage, model_object.GetSelectedAtoms());
}

void RunFixedOffsetLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage,
    const std::vector<AtomObject *> & atom_list)
{
    const auto selected_atom_size{ atom_list.size() };
    auto analysis{ model_object.EditAnalysis() };
    std::vector<LocalGaussianResult> local_results(selected_atom_size);
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
        const auto local_view{ AtomLocalPotentialView::For(atom) };
        LocalPotentialSampleList sample_entries{ local_view.GetSamplingEntries(stage) };
        GaussianModel3D offset_model{ local_view.GetEstimateMDPDE(stage) };
        local_results[i] =
            EstimateLocalGaussian(
                sample_entries,
                local_view.GetAlphaR(stage),
                options,
                offset_model);

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

    for (size_t i = 0; i < selected_atom_size; i++)
    {
        analysis.SetAtomLocalGaussianResult(
            stage, *atom_list[i], std::move(local_results[i]));
    }
}

double TrainAlphaR(
    const LocalPotentialSampleList & sample_entries,
    const FitOptions & options)
{
    std::vector<RHBMMemberDataset> dataset_list{
        rhbm_helper::BuildMemberDataset(
            sample_entries, 0.0, detail::kSignalDistanceMax)
    };
    const auto response_count{
        static_cast<std::size_t>(dataset_list.front().y.size())
    };
    auto training_options{ MakeTrainingOptions(options) };
    if (response_count < 2)
    {
        return training_options.alpha_min;
    }
    training_options.subset_size = std::min(
        training_options.subset_size,
        response_count);
    return rhbm_trainer::CrossValidationAlphaR(dataset_list, training_options).best_alpha;
}

double TrainAlphaG(
    const std::vector<std::vector<GaussianModel3D>> & model_group_list,
    const FitOptions & options)
{
    std::vector<std::vector<RHBMParameterVector>> beta_group_list;
    beta_group_list.reserve(model_group_list.size());
    for (const auto & model_group : model_group_list)
    {
        std::vector<RHBMParameterVector> beta_list;
        beta_list.reserve(model_group.size());
        for (const auto & model : model_group)
        {
            beta_list.emplace_back(
                linearization_service::EncodeGaussianToParameterVector(model));
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
    const detail::PreparedLocalGaussianDesign design{
        sample_entries,
        0.0,
        detail::kSignalDistanceMax
    };
    const auto sample_response_list{ CollectSampleResponses(sample_entries) };
    return design.Estimate(
        sample_response_list,
        alpha_r,
        options.thread_size,
        offset_model);
}

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<GroupGaussianMemberInput> & member_list,
    double alpha_g,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegative(alpha_g, "alpha_g");

    const RHBMExecutionOptions execution_options{
        .thread_size = options.thread_size
    };
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(member_list.size());
    std::vector<RHBMBetaEstimateResult> fit_result_list;
    fit_result_list.reserve(member_list.size());
    std::vector<double> member_offset_list;
    member_offset_list.reserve(member_list.size());
    for (const auto & member : member_list)
    {
        const detail::PreparedLocalGaussianDesign design{
            member.sample_entries,
            0.0,
            detail::kSignalDistanceMax
        };
        auto dataset{
            design.BuildDataset(
                CollectSampleResponses(member.sample_entries),
                member.local_model)
        };
        fit_result_list.emplace_back(
            rhbm_helper::EstimateBetaMDPDE(
                member.alpha_r,
                dataset,
                execution_options));
        dataset_list.emplace_back(std::move(dataset));
        member_offset_list.emplace_back(member.local_model.GetOffset());
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
    model_object.EditAnalysis().EnsureSelectedAtomLocalPotentials();
    RunLocalAlphaTraining(model_object, options, stage, model_object.GetSelectedAtoms());
}

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    FittingStage stage,
    const std::vector<AtomObject *> & atom_list)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto alpha_min{ MakeTrainingOptions(options).alpha_min };

    size_t count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run local alpha training for " +
            std::to_string(atom_list.size()) + " atoms.");
    }
    for (auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        auto alpha_r{ alpha_min };
        if (local_view.HasEnoughSamplingEntriesInRange(
                stage,
                0.0,
                detail::kSignalDistanceMax,
                kMinimumAlphaRTrainingSampleCount))
        {
            alpha_r = TrainAlphaR(
                local_view.GetSamplingEntries(stage),
                options);
        }
        analysis.SetAtomLocalAlphaR(stage, *atom, alpha_r);
        count++;
        if (!options.quiet_mode)
        {
            Logger::ProgressPercent(count, atom_list.size());
        }
    }
}

void RunGroupPotentialFitting(ModelObject & model_object, const FitOptions & options)
{
    if (options.estimator == PotentialEstimator::JOINT_COMPONENTS)
    { detail::RunJointGroupPotentialFitting(model_object, options); return; }
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    analysis.EnsureSelectedAtomLocalPotentials();
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run atom group fitting.");
    }

    auto group_key_list{ analysis_view.CollectAtomGroupKeys() };
    auto group_key_size{ group_key_list.size() };
    size_t key_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t k = 0; k < group_key_size; k++)
    {
        auto group_key{ group_key_list[k] };
        const auto & atom_list{ analysis_view.GetAtomObjectList(group_key) };
        const auto alpha_g{ analysis_view.GetAtomAlphaG(group_key) };
        std::vector<GroupGaussianMemberInput> member_list;
        member_list.reserve(atom_list.size());
        for (const auto & atom : atom_list)
        {
            const auto local_view{ AtomLocalPotentialView::For(*atom) };
            const auto & local_result{ local_view.GetGaussianResult(FittingStage::Second) };
            auto sample_entries{ local_view.GetSamplingEntries(FittingStage::Second) };
            member_list.emplace_back(GroupGaussianMemberInput{
                std::move(sample_entries),
                local_result.alpha_r,
                local_result.mdpde.GetModel()
            });
        }
        const auto result{
            EstimateGroupGaussian(member_list, alpha_g, options)
        };

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            analysis.ApplyAtomGroupGaussianResult(group_key, result);
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
    if (options.estimator != PotentialEstimator::TWO_STAGE)
        throw std::invalid_argument("Joint fitting requires the map-aware workflow.");
    model_object.EditAnalysis().InitializeLocalFittingSeedModels();

    RunLocalAlphaTraining(model_object, options, FittingStage::First);
    RunFixedOffsetLocalFitting(model_object, options, FittingStage::First);

    detail::RunSecondStageIterations(model_object, options);

    if (!options.quiet_mode) Logger::Log(LogLevel::Info, BuildSecondStageSpotSummary(model_object));
    RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options);
}

void RunPotentialFittingWorkflow(MapObject & map, ModelObject & model, const FitOptions & options)
{
    if (options.estimator == PotentialEstimator::TWO_STAGE)
    {
        model.EditAnalysis().InitializeFromSelection();
        RunPotentialSamplingWorkflow(map, model, options.sampling_method, options.thread_size);
        RunPotentialFittingWorkflow(model, options);
        return;
    }
    if (options.sampling_method != SphereSamplingMethod::FibonacciDeterministic)
        throw std::invalid_argument("Joint initialization requires Fibonacci sampling.");
    const auto construction_start = std::chrono::steady_clock::now();
    const auto problem = BuildJointProblem(map, model);
    const auto construction_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - construction_start).count();
    const auto initialization_start = std::chrono::steady_clock::now();
    const auto workset = detail::MakeJointFittingWorkset(model, problem);
    model.EditAnalysis().InitializeFromSelection();
    const auto initialization = detail::RunContributorFirstStage(map, model, workset, options);
    const auto initialization_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - initialization_start).count();
    auto snapshot = [&] {
        auto fit = FitJointComponents(problem, initialization.b);
        fit.costs.construction_seconds = construction_seconds;
        fit.costs.initialization_seconds = initialization_seconds;
        fit.initialization.atoms = initialization.atoms;
        fit.initialization.data_scope = initialization.data_scope;
        return CaptureJointAnalysisResult(fit);
    }();
    data_internal::ApplyJointStageEstimates(model, snapshot, boost::uuids::to_string(boost::uuids::random_generator()()));
    if (!options.quiet_mode) Logger::Log(LogLevel::Info, BuildSecondStageSpotSummary(model));
    for (auto & [id, peeling] : detail::BuildPostFitPeelingSamples(map, model, problem, problem.Input().selection_domain->target_indices))
        model.EditAnalysis().SetAtomPostFitPeeling(*model.FindAtomPtr(id), std::move(peeling));
    for (auto & [id, uncertainty] : detail::ComputeJointUncertainty(problem, snapshot, problem.Input().selection_domain->target_indices))
    {
        auto & atom = *model.FindAtomPtr(id);
        auto stage = AtomLocalPotentialView::For(atom).GetStageEstimate(FittingStage::Second);
        stage.uncertainty = std::move(uncertainty);
        model.EditAnalysis().SetAtomStageEstimate(FittingStage::Second, atom, stage);
        model.EditAnalysis().SetAtomGroupEvidence(atom, detail::BuildJointParameterEvidence(stage));
    }
    model.EditAnalysis().SetJointResult(std::move(snapshot));
    RunGroupPotentialFitting(model, options);
}

} // namespace rhbm_gem::core
