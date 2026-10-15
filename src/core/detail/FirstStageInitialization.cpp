#include "FirstStageInitialization.hpp"
#include "FittingWorkset.hpp"
#include "gaussian_fit/FittingRanges.hpp"
#include "gaussian_fit/PreparedLocalGaussianFit.hpp"
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include "utils/hrl/MDPDEEndpointRefinement.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace rhbm_gem::core::detail {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
FirstStageObserver & FirstStageObserverForTesting()
{
    static thread_local FirstStageObserver observer;
    return observer;
}
#endif
void ApplyJointSeedFallback(JointInitialization & initialization)
{
    std::vector<double> donors;
    for(const auto & atom:initialization.atoms)
        if(atom.target && atom.reason=="valid-width" && atom.original_b &&
            std::isfinite(*atom.original_b) && *atom.original_b>0)
            donors.push_back(*atom.original_b);
    std::sort(donors.begin(),donors.end());
    const double median=donors.empty() ? 0 : donors.size()%2 ? donors[donors.size()/2] :
        donors[donors.size()/2-1]+(donors[donors.size()/2]-donors[donors.size()/2-1])/2;
    initialization.valid=true;
    for(std::size_t i=0;i<initialization.atoms.size();++i)
    {
        auto & atom=initialization.atoms[i];
        if(atom.target && !(std::isfinite(initialization.b[i]) && initialization.b[i]>0) && !donors.empty())
        {initialization.b[i]=median; atom.seed_source="median-fallback"; atom.donor_count=donors.size();}
        const bool valid=std::isfinite(initialization.b[i]) && initialization.b[i]>0;
        if(valid) atom.used_b=initialization.b[i];
        if(atom.target) initialization.valid &= valid;
        else
        {
            std::string failure;
            if(!MakeFixedHaloParameter(atom,i,atom.id,failure))
            {
                atom.reason=failure;
                initialization.valid=false;
                initialization.reason="fixed-halo-initialization-unavailable:"+failure;
            }
        }
    }
    if(initialization.reason.empty())
        initialization.reason=initialization.valid ? "valid-widths" : "invalid-widths";
}

std::optional<JointFixedHaloParameter> MakeFixedHaloParameter(
    const JointInitializationAtom & atom, std::size_t contributor_index,
    std::string_view expected_atom_id, std::string & failure_reason)
{
    failure_reason.clear();
    if (expected_atom_id.empty() || atom.id.empty())
        failure_reason = "missing-atom-identity";
    else if (atom.id != expected_atom_id)
        failure_reason = "contributor-identity-mismatch";
    else if (atom.target)
        failure_reason = "target-is-not-fixed-halo";
    else if (atom.seed_source == "median-fallback")
        failure_reason = "median-width-fallback";
    else if (!IsQualifiedFirstStageEstimate(atom.native_status, atom.first_stage_fit_qualified,
        atom.first_stage_qualification, atom.first_stage_recovery))
        failure_reason = "first-stage-fitting-failure";
    else if (!std::isfinite(atom.mdpde[0]) || !std::isfinite(atom.mdpde[1]))
        failure_reason = "non-finite-parameters";
    else if (atom.mdpde[0] < 0.0)
        failure_reason = "invalid-amplitude";
    else if (!(atom.mdpde[1] > 0.0))
        failure_reason = "invalid-width";
    if (!failure_reason.empty()) return std::nullopt;
    JointFixedHaloParameter parameter;
    parameter.contributor_index = contributor_index;
    parameter.atom_id = atom.id;
    parameter.amplitude = atom.mdpde[0];
    parameter.width = atom.mdpde[1];
    parameter.charge = 0.0;
    parameter.first_stage_sample_count = atom.sample_count;
    parameter.first_stage_status = atom.native_status;
    parameter.first_stage_fit_qualified = atom.first_stage_fit_qualified;
    parameter.first_stage_qualification = FirstStageQualification(atom.native_status,
        atom.first_stage_fit_qualified, atom.first_stage_qualification, atom.first_stage_recovery);
    parameter.first_stage_recovery = atom.first_stage_recovery;
    return parameter;
}

RHBMBetaEstimateResult RecoverFailedFirstStageEstimate(const RHBMMemberDataset & dataset, double alpha,
    const RHBMExecutionOptions & native_options, const RHBMBetaEstimateResult & native)
{
    if (native.status == RHBMEstimationStatus::SUCCESS) return native;
    RHBMFailedOnlyRecoveryDiagnostics recovery;
    const auto refinement_start = std::chrono::steady_clock::now();
    auto refined = mdpde_detail::ApplyFailedOnlyRefinement(dataset, alpha, native_options, native);
    recovery.endpoint_refinement_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - refinement_start).count();
    recovery.endpoint_refinement = refined.refinement;
    auto qualified = refined;
    if (refined.Qualification() == RHBMSolveQualification::RefinedSuccess)
    {
        recovery.accepted = true;
        recovery.reason = "endpoint-refinement-accepted";
    }
    else if (native.status == RHBMEstimationStatus::MAX_ITERATIONS_REACHED)
    {
        auto retry_options = native_options;
        retry_options.max_iterations = 300;
        const auto retry_start = std::chrono::steady_clock::now();
        const auto retried = rhbm_helper::EstimateBetaMDPDE(alpha, dataset, retry_options);
        RHBMNativeRetryDiagnostics retry;
        retry.max_iterations = retry_options.max_iterations;
        retry.native_status = static_cast<int>(retried.status);
        retry.iterations = retried.diagnostics.iterations;
        retry.squared_beta_change = retried.diagnostics.squared_beta_change;
        retry.relative_variance_change = retried.diagnostics.relative_variance_change;
        const auto equations = mdpde_detail::EvaluateMDPDEEquations(dataset, alpha,
            retried.beta_mdpde, retried.sigma_square, retry_options.data_weight_min);
        if (equations.valid && equations.scaled.allFinite())
            retry.equation_residual = equations.scaled.lpNorm<Eigen::Infinity>();
        retry.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - retry_start).count();
        recovery.bounded_native_retry = retry;
        if (retried.status == RHBMEstimationStatus::SUCCESS)
        {
            qualified = retried;
            qualified.status = native.status;
            qualified.diagnostics = native.diagnostics;
            qualified.refinement = refined.refinement;
            recovery.accepted = true;
            recovery.reason = "bounded-native-retry-accepted";
        }
        else
        {
            qualified = native;
            qualified.refinement = refined.refinement;
            recovery.reason = "bounded-native-retry-unqualified";
        }
    }
    else recovery.reason = refined.refinement ? refined.refinement->reason : "endpoint-refinement-unavailable";

    qualified.failed_only_recovery = std::move(recovery);
    return qualified;
}

LocalGaussianResult FitFirstStageAtom(const AtomObject & atom, const FitOptions & options,
    bool enable_fixed_halo_recovery)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if (FirstStageObserverForTesting()) FirstStageObserverForTesting()(atom.GetSerialID(), "first");
#endif
    const auto view = AtomLocalPotentialView::For(atom);
    const auto samples = view.GetSamplingEntries(FittingStage::First);
    auto result = EstimateLocalGaussian(samples, view.GetAlphaR(FittingStage::First),
        options, view.GetEstimateMDPDE(FittingStage::First));
    if (!enable_fixed_halo_recovery || !result.fit_result ||
        result.fit_result->status == RHBMEstimationStatus::SUCCESS)
        return result;

    std::vector<double> responses;
    responses.reserve(samples.size());
    for (const auto & sample : samples) responses.push_back(sample.response);
    const PreparedLocalGaussianDesign design{samples, 0.0, kSignalDistanceMax};
    const auto dataset = design.BuildDataset(responses, view.GetEstimateMDPDE(FittingStage::First));
    RHBMExecutionOptions native_options{.thread_size = options.thread_size};
    auto qualified = RecoverFailedFirstStageEstimate(dataset, view.GetAlphaR(FittingStage::First),
        native_options, *result.fit_result);
    result.mdpde = GaussianModel3DWithUncertainty{
        linearization_service::DecodeParameterVector(qualified.beta_mdpde), GaussianModel3DUncertainty{}};
    result.fit_result = std::move(qualified);
    result.diagnostics = LocalFitDiagnostics{result.fit_result->status, result.fit_result->sigma_square,
        result.fit_result->diagnostics, result.fit_result->refinement, result.fit_result->failed_only_recovery};
    return result;
}

void RunBatchFirstStageFromPreparedSamples(
    ModelObject & model, const std::vector<AtomObject *> & atoms, const FitOptions & options)
{
    RunLocalAlphaTraining(model, options, FittingStage::First, atoms);
    RunFixedOffsetLocalFitting(model, options, FittingStage::First, atoms);
}

JointInitialization RunJointFirstStageInitializationFromPreparedSamples(
    ModelObject & model, const FittingWorkset & workset, const FitOptions & options)
{
    if (workset.target_mask.size() != workset.contributors.size())
        throw std::invalid_argument("First-stage workset role count mismatch.");
    JointInitialization initialization;
    initialization.data_scope = "contributor-local-sampling-may-read-outside-target-domain";
    auto editor = model.EditAnalysis();
    std::size_t completed_count = 0;
    for (std::size_t index = 0; index < workset.contributors.size(); ++index)
    {
        auto * atom = workset.contributors[index];
        JointInitializationAtom record;
        record.id = std::to_string(atom->GetSerialID());
        record.target = workset.target_mask[index];
        record.ols.fill(std::numeric_limits<double>::quiet_NaN());
        record.mdpde = record.ols;
        record.alpha = std::numeric_limits<double>::quiet_NaN();
        double width = std::numeric_limits<double>::quiet_NaN();
        try
        {
            LocalGaussianResult seed;
            seed.ols = seed.mdpde = GaussianModel3DWithUncertainty{GaussianModel3D{0, 1, 0}, {}};
            editor.SetAtomLocalGaussianResult(FittingStage::First, *atom, seed);
            editor.SetAtomStageEstimate(FittingStage::First, *atom, LocalStageEstimate{});
            const auto view = AtomLocalPotentialView::For(*atom);
            record.sample_count = view.GetSamplingEntries(FittingStage::First).size();
            TrainLocalAlphaForAtom(editor, options, FittingStage::First, *atom);
            record.alpha = view.GetAlphaR(FittingStage::First);
            editor.SetAtomLocalGaussianResult(FittingStage::First, *atom,
                FitFirstStageAtom(*atom, options, !workset.target_mask[index]));
            auto first = view.GetStageEstimate(FittingStage::First);
            first.source.role = workset.target_mask[index] ? FittingRole::Target : FittingRole::Halo;
            editor.SetAtomStageEstimate(FittingStage::First, *atom, first);
            const auto & local = view.GetGaussianResult(FittingStage::First);
            const auto ols = local.ols.GetModel().ToVector(), mdpde = local.mdpde.GetModel().ToVector();
            for (std::size_t k = 0; k < 3; ++k) { record.ols[k] = ols(static_cast<Eigen::Index>(k)); record.mdpde[k] = mdpde(static_cast<Eigen::Index>(k)); }
            if (local.fit_result)
            {
                record.native_status = static_cast<int>(local.fit_result->status);
                record.first_stage_qualification = local.fit_result->Qualification();
                record.first_stage_recovery = local.fit_result->failed_only_recovery;
                record.first_stage_fit_qualified = record.first_stage_qualification !=
                    RHBMSolveQualification::Unqualified;
            }
            width = local.mdpde.GetModel().GetWidth();
            if (!record.target && record.first_stage_qualification == RHBMSolveQualification::Unqualified)
            {
                if (record.first_stage_recovery)
                    record.reason = "first-stage-recovery-unqualified:" + record.first_stage_recovery->reason;
                else if (record.native_status)
                    record.reason = "first-stage-native-unqualified:" + std::to_string(*record.native_status);
                else
                    record.reason = "first-stage-result-unavailable";
            }
            else record.reason = std::isfinite(width) && width > 0 ? "valid-width" : "invalid-width";
        }
        catch (const std::exception & error)
        {
            LocalStageEstimate missing;
            missing.reason = std::string("initialization-exception: ") + error.what();
            editor.SetAtomStageEstimate(FittingStage::First, *atom, missing);
            record.reason = missing.reason;
        }
        if (std::isfinite(width)) record.original_b = width;
        record.seed_source = std::isfinite(width) && width > 0 ? "fitted" : "unavailable";
        initialization.b.push_back(width);
        initialization.atoms.push_back(std::move(record));
        ++completed_count;
        if (!options.quiet_mode)
        {
            Logger::ProgressPercent(completed_count, workset.contributors.size(), 50, " Joint initialization");
        }
    }
    ApplyJointSeedFallback(initialization);
    return initialization;
}

JointInitialization PrepareJointWidthInitialization(
    MapObject & map, ModelObject & model, const FittingWorkset & workset, const FitOptions & options)
{
    if (!options.quiet_mode)
    {
        Logger::FinishProgressLine();
        Logger::Log(LogLevel::Info, "[Joint] Sampling " +
            std::to_string(workset.contributors.size()) +
            " FullABC contributors");
    }
    const auto atoms = CollectFirstStageAtoms(workset);
    RunPotentialSamplingWorkflow(map, model, atoms, SphereSamplingMethod::FibonacciDeterministic,
        1, " Joint sampling", options.quiet_mode);
    if (!options.quiet_mode)
    {
        Logger::FinishProgressLine();
        Logger::Log(LogLevel::Info, "[Joint] Initializing " +
            std::to_string(workset.contributors.size()) + " contributors");
    }
    return RunJointFirstStageInitializationFromPreparedSamples(model, workset, options);
}
}
