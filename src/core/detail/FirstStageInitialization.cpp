#include "FirstStageInitialization.hpp"
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
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
FittingWorkset MakeJointFittingWorkset(ModelObject & model, const JointProblem & problem)
{
    const auto & input = problem.Input();
    if (!input.selection_domain) throw std::invalid_argument("Joint workset requires recorded selection.");
    FittingWorkset workset;
    for (std::size_t i = 0; i < input.atom_ids.size(); ++i)
    {
        auto * atom = model.FindAtomPtr(std::stoi(input.atom_ids[i]));
        if (!atom || std::to_string(atom->GetSerialID()) != input.atom_ids[i])
            throw std::invalid_argument("Unknown joint contributor identity.");
        workset.contributors.push_back(atom);
        const auto & targets = input.selection_domain->target_indices;
        workset.target_mask.push_back(std::binary_search(targets.begin(), targets.end(), i));
        const auto & full=problem.ParameterLayout().full_atoms;
        workset.full_parameter_mask.push_back(std::binary_search(full.begin(),full.end(),i));
    }
    return workset;
}

void ApplyJointSeedFallback(JointInitialization & initialization)
{
    std::vector<double> donors;
    for(const auto & atom:initialization.atoms)
        if(atom.reason=="valid-width" && atom.original_b && std::isfinite(*atom.original_b) && *atom.original_b>0)
            donors.push_back(*atom.original_b);
    std::sort(donors.begin(),donors.end());
    const double median=donors.empty() ? 0 : donors.size()%2 ? donors[donors.size()/2] :
        donors[donors.size()/2-1]+(donors[donors.size()/2]-donors[donors.size()/2-1])/2;
    initialization.valid=true;
    for(std::size_t i=0;i<initialization.atoms.size();++i)
    {
        auto & atom=initialization.atoms[i];
        if(atom.reason=="not-required-observable-contribution") continue;
        if(!(std::isfinite(initialization.b[i]) && initialization.b[i]>0) && !donors.empty())
        {initialization.b[i]=median; atom.seed_source="median-fallback"; atom.donor_count=donors.size();}
        const bool valid=std::isfinite(initialization.b[i]) && initialization.b[i]>0;
        if(valid) atom.used_b=initialization.b[i];
        initialization.valid &= valid;
    }
    initialization.reason=initialization.valid ? "valid-widths" : "invalid-widths";
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

JointInitialization RunFirstStage(ModelObject & model, const FittingWorkset & workset,
    const FitOptions & options, FirstStageMode mode, MapObject * sampling_map)
{
    if (workset.target_mask.size() != workset.contributors.size())
        throw std::invalid_argument("First-stage workset role count mismatch.");
    if (mode == FirstStageMode::ExistingSamplesBatch)
    {
        RunLocalAlphaTraining(model, options, FittingStage::First, workset.contributors);
        RunFixedOffsetLocalFitting(model, options, FittingStage::First, workset.contributors);
        return {};
    }
    if (!sampling_map) throw std::invalid_argument("Contributor sampling requires a map.");
    auto & map = *sampling_map;
    JointInitialization initialization;
    initialization.data_scope = "contributor-local-sampling-may-read-outside-target-domain";
    auto editor = model.EditAnalysis();
    for (std::size_t index = 0; index < workset.contributors.size(); ++index)
    {
        auto * atom = workset.contributors[index];
        JointInitializationAtom record;
        record.id = std::to_string(atom->GetSerialID());
        record.ols.fill(std::numeric_limits<double>::quiet_NaN());
        record.mdpde = record.ols;
        record.alpha = std::numeric_limits<double>::quiet_NaN();
        double width = std::numeric_limits<double>::quiet_NaN();
        if(!workset.full_parameter_mask.empty() && !workset.full_parameter_mask.at(index))
        {
            record.reason="not-required-observable-contribution"; record.seed_source="not-required";
            initialization.b.push_back(width); initialization.atoms.push_back(std::move(record)); continue;
        }
        try
        {
            LocalGaussianResult seed;
            seed.ols = seed.mdpde = GaussianModel3DWithUncertainty{GaussianModel3D{0, 1, 0}, {}};
            editor.SetAtomLocalGaussianResult(FittingStage::First, *atom, seed);
            editor.SetAtomStageEstimate(FittingStage::First, *atom, LocalStageEstimate{});
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if (FirstStageObserverForTesting()) FirstStageObserverForTesting()(atom->GetSerialID(), "raw");
#endif
            editor.SetAtomLocalRawSamplingEntries(*atom, SampleAtomMapValues(map, *atom, options.sampling_method));
            const auto view = AtomLocalPotentialView::For(*atom);
            record.sample_count = view.GetSamplingEntries(FittingStage::First).size();
            const std::vector<AtomObject *> atoms{atom};
            RunLocalAlphaTraining(model, options, FittingStage::First, atoms);
            record.alpha = view.GetAlphaR(FittingStage::First);
            editor.SetAtomLocalGaussianResult(FittingStage::First, *atom, FitFirstStageAtom(*atom, options));
            auto first = view.GetStageEstimate(FittingStage::First);
            first.source.role = workset.target_mask[index] ? FittingRole::Target : FittingRole::Halo;
            editor.SetAtomStageEstimate(FittingStage::First, *atom, first);
            const auto & local = view.GetGaussianResult(FittingStage::First);
            const auto ols = local.ols.GetModel().ToVector(), mdpde = local.mdpde.GetModel().ToVector();
            for (std::size_t k = 0; k < 3; ++k) { record.ols[k] = ols(static_cast<Eigen::Index>(k)); record.mdpde[k] = mdpde(static_cast<Eigen::Index>(k)); }
            if (local.fit_result) record.native_status = static_cast<int>(local.fit_result->status);
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
        initialization.b.push_back(width);
        if(std::isfinite(width)) record.original_b=width;
        record.seed_source=std::isfinite(width) && width>0 ? "fitted" : "unavailable";
        initialization.atoms.push_back(std::move(record));
    }
    ApplyJointSeedFallback(initialization);
    return initialization;
}
}
