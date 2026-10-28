#include "FirstStageInitialization.hpp"
#include "FittingWorkset.hpp"
#include <rhbm_gem/core/MapSampler.hpp>
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
    }
    if(initialization.reason.empty())
        initialization.reason=initialization.valid ? "valid-widths" : "invalid-widths";
}

LocalGaussianResult FitFirstStageAtom(const AtomObject & atom, const FitOptions & options)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if (FirstStageObserverForTesting()) FirstStageObserverForTesting()(atom.GetSerialID(), "first");
#endif
    const auto view = AtomLocalPotentialView::For(atom);
    const auto samples = view.GetSamplingEntries(FittingStage::First);
    return EstimateLocalGaussian(samples, view.GetAlphaR(FittingStage::First),
        options, view.GetEstimateMDPDE(FittingStage::First));
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
    std::optional<ScopeTimer> timer;
    if (!options.quiet_mode) timer.emplace("Joint FirstStage initialization");
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    const auto worker_count = std::clamp(options.joint_initialization_thread_size, 1, 4);
    FitOptions atom_options = options;
    atom_options.thread_size = 1;
    initialization.atoms.resize(workset.contributors.size());
    initialization.b.assign(workset.contributors.size(), nan);
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
        editor.SetAtomLocalGaussianResult(FittingStage::First, *atom, seed);
        editor.SetAtomStageEstimate(FittingStage::First, *atom, LocalStageEstimate{});
    }

    const auto workspace_count = static_cast<std::size_t>(worker_count);
    std::vector<LocalGaussianResult> estimates(workspace_count);
    std::vector<double> alphas(workspace_count, nan);
    std::vector<std::size_t> sample_counts(workspace_count);
    std::vector<unsigned char> alpha_ready(workspace_count);
    std::vector<unsigned char> fit_attempted(workspace_count);
    std::vector<std::exception_ptr> failures(workspace_count);
    std::size_t completed_count = 0;
    for (std::size_t begin = 0; begin < workset.contributors.size(); begin += workspace_count)
    {
        const auto count = std::min(workspace_count, workset.contributors.size() - begin);
        std::fill(alphas.begin(), alphas.end(), nan);
        std::fill(sample_counts.begin(), sample_counts.end(), 0);
        std::fill(alpha_ready.begin(), alpha_ready.end(), 0);
        std::fill(fit_attempted.begin(), fit_attempted.end(), 0);
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
                const auto & samples = view.GetSamplingEntries(FittingStage::First);
                sample_counts[slot] = samples.size();
                alphas[slot] = EstimateLocalAlphaForAtom(*atom, atom_options, FittingStage::First);
                alpha_ready[slot] = 1;
                fit_attempted[slot] = 1;
                estimates[slot] = EstimateLocalGaussian(samples, alphas[slot], atom_options,
                    view.GetEstimateMDPDE(FittingStage::First));
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
            double width = nan;
            try
            {
                if (alpha_ready[slot])
                {
                    editor.SetAtomLocalAlphaR(FittingStage::First, *atom, alphas[slot]);
                    record.alpha = alphas[slot];
                }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
                if (fit_attempted[slot] && FirstStageObserverForTesting())
                    FirstStageObserverForTesting()(atom->GetSerialID(), "first");
#endif
                if (failures[slot]) std::rethrow_exception(failures[slot]);
                editor.SetAtomLocalGaussianResult(FittingStage::First, *atom, std::move(estimates[slot]));
                const auto view = AtomLocalPotentialView::For(*atom);
                auto first = view.GetStageEstimate(FittingStage::First);
                first.source.role = record.target ? FittingRole::Target : FittingRole::Halo;
                editor.SetAtomStageEstimate(FittingStage::First, *atom, first);
                const auto & local = view.GetGaussianResult(FittingStage::First);
                const auto ols = local.ols.GetModel().ToVector(), mdpde = local.mdpde.GetModel().ToVector();
                for (std::size_t k = 0; k < 3; ++k)
                {
                    record.ols[k] = ols(static_cast<Eigen::Index>(k));
                    record.mdpde[k] = mdpde(static_cast<Eigen::Index>(k));
                }
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
            initialization.b[index] = width;
            ++completed_count;
            if (!options.quiet_mode)
            {
                Logger::ProgressPercent(completed_count, workset.contributors.size(), 50, " Joint initialization");
            }
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
            " Joint contributors");
    }
    const auto atoms = CollectFirstStageAtoms(workset);
    RunPotentialSamplingWorkflow(map, model, atoms, options.joint_initialization_thread_size,
        " Joint sampling", options.quiet_mode);
    if (!options.quiet_mode)
    {
        Logger::FinishProgressLine();
        Logger::Log(LogLevel::Info, "[Joint] Initializing " +
            std::to_string(workset.contributors.size()) + " contributors");
    }
    return RunJointFirstStageInitializationFromPreparedSamples(model, workset, options);
}
}
