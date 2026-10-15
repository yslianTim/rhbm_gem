#include "FirstStageInitialization.hpp"
#include "FittingWorkset.hpp"
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

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
    else if (!atom.native_status || !atom.first_stage_fit_qualified)
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
    return parameter;
}

LocalGaussianResult FitFirstStageAtom(const AtomObject & atom, const FitOptions & options)
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if (FirstStageObserverForTesting()) FirstStageObserverForTesting()(atom.GetSerialID(), "first");
#endif
    const auto view = AtomLocalPotentialView::For(atom);
    return EstimateLocalGaussian(view.GetSamplingEntries(FittingStage::First), view.GetAlphaR(FittingStage::First),
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
            editor.SetAtomLocalGaussianResult(FittingStage::First, *atom, FitFirstStageAtom(*atom, options));
            auto first = view.GetStageEstimate(FittingStage::First);
            first.source.role = workset.target_mask[index] ? FittingRole::Target : FittingRole::Halo;
            editor.SetAtomStageEstimate(FittingStage::First, *atom, first);
            const auto & local = view.GetGaussianResult(FittingStage::First);
            const auto ols = local.ols.GetModel().ToVector(), mdpde = local.mdpde.GetModel().ToVector();
            for (std::size_t k = 0; k < 3; ++k) { record.ols[k] = ols(static_cast<Eigen::Index>(k)); record.mdpde[k] = mdpde(static_cast<Eigen::Index>(k)); }
            if (local.fit_result)
            {
                record.native_status = static_cast<int>(local.fit_result->status);
                record.first_stage_fit_qualified =
                    local.fit_result->status == RHBMEstimationStatus::SUCCESS &&
                    local.fit_result->Qualification() != RHBMSolveQualification::Unqualified;
            }
            width = local.mdpde.GetModel().GetWidth();
            record.reason = std::isfinite(width) && width > 0 ? "valid-width" : "invalid-width";
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
