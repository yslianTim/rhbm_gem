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
    }
    return workset;
}

JointInitialization RunContributorFirstStage(MapObject & map, ModelObject & model,
    const FittingWorkset & workset, const FitOptions & options)
{
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
            if (FirstStageObserverForTesting()) FirstStageObserverForTesting()(atom->GetSerialID(), "first");
#endif
            editor.SetAtomLocalGaussianResult(FittingStage::First, *atom,
                EstimateLocalGaussian(view.GetSamplingEntries(FittingStage::First), record.alpha,
                    options, seed.mdpde.GetModel()));
            auto first = view.GetStageEstimate(FittingStage::First);
            first.source.role = workset.target_mask[index] ? FittingRole::Target : FittingRole::Halo;
            editor.SetAtomStageEstimate(FittingStage::First, *atom, first);
            const auto & local = view.GetGaussianResult(FittingStage::First);
            const auto ols = local.ols.GetModel().ToVector(), mdpde = local.mdpde.GetModel().ToVector();
            for (std::size_t k = 0; k < 3; ++k) { record.ols[k] = ols(k); record.mdpde[k] = mdpde(k); }
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
        initialization.atoms.push_back(std::move(record));
    }
    initialization.valid = std::all_of(initialization.b.begin(), initialization.b.end(),
        [](double b) { return std::isfinite(b) && b > 0; });
    initialization.reason = initialization.valid ? "valid-widths" : "invalid-widths";
    return initialization;
}
}
