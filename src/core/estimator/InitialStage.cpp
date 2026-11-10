#include "InitialStage.hpp"

#include "core/sampling/MapSampler.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace rhbm_gem::core::detail {
InitialStageWorkset BuildInitialStageWorkset(ModelObject & model, const JointProblem & problem)
{
    const auto & input = problem.Input();
    if (!input.selection_domain) throw std::invalid_argument("Joint workset requires recorded selection.");
    InitialStageWorkset workset;
    for (std::size_t i = 0; i < input.atom_ids.size(); ++i)
    {
        auto * atom = model.FindAtomPtr(std::stoi(input.atom_ids[i]));
        if (!atom || std::to_string(atom->GetSerialID()) != input.atom_ids[i])
            throw std::invalid_argument("Unknown joint contributor identity.");
        workset.contributors.push_back(atom);
        const auto & targets = input.selection_domain->target_indices;
        workset.target_mask.push_back(std::binary_search(targets.begin(), targets.end(), i));
    }
    return workset;
}

void ApplyJointSeedFallback(JointInitialization & initialization)
{
    std::vector<double> donors;
    for(const auto & atom:initialization.atoms)
        if(atom.target && atom.reason=="valid-width" && atom.original_tau &&
            std::isfinite(*atom.original_tau) && *atom.original_tau>0)
            donors.push_back(*atom.original_tau);
    std::sort(donors.begin(),donors.end());
    const double median=donors.empty() ? 0 : donors.size()%2 ? donors[donors.size()/2] :
        donors[donors.size()/2-1]+(donors[donors.size()/2]-donors[donors.size()/2-1])/2;
    initialization.valid=true;
    for(std::size_t i=0;i<initialization.atoms.size();++i)
    {
        auto & atom=initialization.atoms[i];
        if(atom.target && !(std::isfinite(initialization.tau[i]) && initialization.tau[i]>0) && !donors.empty())
        {initialization.tau[i]=median; atom.seed_source="median-fallback"; atom.donor_count=donors.size();}
        const bool valid=std::isfinite(initialization.tau[i]) && initialization.tau[i]>0;
        if(valid) atom.used_tau=initialization.tau[i];
        if(atom.target) initialization.valid &= valid;
    }
    if(initialization.reason.empty())
        initialization.reason=initialization.valid ? "valid-widths" : "invalid-widths";
}

LocalGaussianResult FitInitialStageAtom(const AtomObject & atom, const FitOptions & options)
{
    const auto view = AtomLocalPotentialView::For(atom);
    const auto samples = view.GetRawSamplingEntries();
    return EstimateLocalGaussian(samples, view.GetInitialStageAlphaR(),
        options, view.GetInitialStageMDPDE());
}

JointInitialization FitInitialStageFromPreparedSamples(
    ModelObject & model, const InitialStageWorkset & workset, const FitOptions & options)
{
    if (workset.target_mask.size() != workset.contributors.size())
        throw std::invalid_argument("InitialStage workset role count mismatch.");
    JointInitialization initialization;
    initialization.data_scope = "contributor-local-sampling-may-read-outside-target-domain";
    auto editor = model.EditAnalysis();
    std::optional<ScopeTimer> timer;
    if (!options.quiet_mode) timer.emplace("[InitialStage] Initialization");
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    const auto worker_count = std::clamp(options.joint_initialization_thread_size, 1, 4);
    FitOptions atom_options = options;
    atom_options.thread_size = 1;
    initialization.atoms.resize(workset.contributors.size());
    initialization.tau.assign(workset.contributors.size(), nan);
    for (std::size_t index = 0; index < workset.contributors.size(); ++index)
    {
        auto * atom = workset.contributors[index];
        auto & record = initialization.atoms[index];
        record.id = std::to_string(atom->GetSerialID());
        record.target = workset.target_mask[index];
        record.ols.fill(nan);
        record.mdpde = record.ols;
        record.alpha = nan;
        LocalGaussianResult seed;
        seed.ols = seed.mdpde = GaussianModel3DWithUncertainty{GaussianModel3D{0, 1, 0}, {}};
        editor.SetAtomInitialStageGaussianResult(*atom, seed);
        editor.SetAtomInitialStageEstimate(*atom, LocalStageEstimate{});
    }

    const auto workspace_count = static_cast<std::size_t>(worker_count);
    std::vector<LocalGaussianResult> estimates(workspace_count);
    std::vector<double> alphas(workspace_count, nan);
    std::vector<std::size_t> sample_counts(workspace_count);
    std::vector<unsigned char> alpha_ready(workspace_count);
    std::vector<std::exception_ptr> failures(workspace_count);
    std::size_t completed_count = 0;
    for (std::size_t begin = 0; begin < workset.contributors.size(); begin += workspace_count)
    {
        const auto count = std::min(workspace_count, workset.contributors.size() - begin);
        std::fill(alphas.begin(), alphas.end(), nan);
        std::fill(sample_counts.begin(), sample_counts.end(), 0);
        std::fill(alpha_ready.begin(), alpha_ready.end(), 0);
        std::fill(failures.begin(), failures.end(), std::exception_ptr{});
#ifdef USE_OPENMP
        #pragma omp parallel for num_threads(worker_count)
#endif
        for (std::ptrdiff_t offset = 0; offset < static_cast<std::ptrdiff_t>(count); ++offset)
        {
            const auto slot = static_cast<std::size_t>(offset);
            auto * atom = workset.contributors[begin + slot];
            try
            {
                const auto view = AtomLocalPotentialView::For(*atom);
                const auto samples = view.GetRawSamplingEntries();
                sample_counts[slot] = samples.size();
                alphas[slot] = EstimateInitialStageAlpha(samples, atom_options);
                alpha_ready[slot] = 1;
                estimates[slot] = EstimateLocalGaussian(samples, alphas[slot], atom_options,
                    view.GetInitialStageMDPDE());
            }
            catch (...)
            {
                failures[slot] = std::current_exception();
            }
        }
        for (std::size_t slot = 0; slot < count; ++slot)
        {
            const auto index = begin + slot;
            auto * atom = workset.contributors[index];
            auto & record = initialization.atoms[index];
            record.sample_count = sample_counts[slot];
            double tau = nan;
            try
            {
                if (alpha_ready[slot])
                {
                    editor.SetAtomInitialStageAlphaR(*atom, alphas[slot]);
                    record.alpha = alphas[slot];
                }
                if (failures[slot]) std::rethrow_exception(failures[slot]);
                editor.SetAtomInitialStageGaussianResult(*atom, std::move(estimates[slot]));
                const auto view = AtomLocalPotentialView::For(*atom);
                auto estimate = view.GetInitialStageEstimate();
                estimate.source.role = record.target ? FittingRole::Target : FittingRole::Halo;
                editor.SetAtomInitialStageEstimate(*atom, estimate);
                const auto & local = view.GetInitialStageGaussianResult();
                const auto ols = local.ols.GetModel().ToVector(), mdpde = local.mdpde.GetModel().ToVector();
                for (std::size_t k = 0; k < 3; ++k)
                {
                    record.ols[k] = ols(static_cast<Eigen::Index>(k));
                    record.mdpde[k] = mdpde(static_cast<Eigen::Index>(k));
                }
                if (local.fit_result)
                {
                    record.native_status = static_cast<int>(local.fit_result->status);
                    record.initial_stage_qualification = local.fit_result->Qualification();
                    record.initial_stage_recovery = local.fit_result->failed_only_recovery;
                    record.initial_stage_fit_qualified = record.initial_stage_qualification !=
                        RHBMSolveQualification::Unqualified;
                }
                tau = local.mdpde.GetModel().GetWidth();
                if (!record.target && record.initial_stage_qualification == RHBMSolveQualification::Unqualified)
                {
                    if (record.initial_stage_recovery)
                        record.reason = "initial-stage-recovery-unqualified:" + record.initial_stage_recovery->reason;
                    else if (record.native_status)
                        record.reason = "initial-stage-native-unqualified:" + std::to_string(*record.native_status);
                    else
                        record.reason = "initial-stage-result-unavailable";
                }
                else record.reason = std::isfinite(tau) && tau > 0 ? "valid-width" : "invalid-width";
            }
            catch (const std::exception & error)
            {
                LocalStageEstimate missing;
                missing.reason = std::string("initialization-exception: ") + error.what();
                editor.SetAtomInitialStageEstimate(*atom, missing);
                record.reason = missing.reason;
            }
            if (std::isfinite(tau)) record.original_tau = tau;
            record.seed_source = std::isfinite(tau) && tau > 0 ? "fitted" : "unavailable";
            initialization.tau[index] = tau;
            ++completed_count;
            if (!options.quiet_mode)
            {
                Logger::ProgressPercent(completed_count, workset.contributors.size(), 50, " [InitialStage] initialization");
            }
        }
    }
    ApplyJointSeedFallback(initialization);
    return initialization;
}

JointInitialization RunInitialStage(
    MapObject & map, ModelObject & model, const InitialStageWorkset & workset, const FitOptions & options)
{
    if (!options.quiet_mode)
    {
        const auto sampling_name{ options.initial_stage_sampling_method == SphereSamplingMethod::SpiralDeterministic
            ? "SpiralShells" : "FibonacciShells" };
        Logger::FinishProgressLine();
        Logger::Log(LogLevel::Info, "[InitialStage] Sampling " +
            std::to_string(workset.contributors.size()) +
            " Joint contributors | sampling=" + sampling_name + " | interpolation=Tricubic");
    }
    RunPotentialSamplingWorkflow(map, model, workset.contributors, std::clamp(options.joint_initialization_thread_size, 1, 4),
        " [InitialStage] sampling", options.quiet_mode, options.initial_stage_sampling_method);
    if (!options.quiet_mode)
    {
        Logger::FinishProgressLine();
        Logger::Log(LogLevel::Info, "[InitialStage] Initializing " +
            std::to_string(workset.contributors.size()) + " contributors");
    }
    return FitInitialStageFromPreparedSamples(model, workset, options);
}
}
