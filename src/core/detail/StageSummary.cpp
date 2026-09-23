#include "StageSummary.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <algorithm>
#include <array>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

namespace rhbm_gem::core {
namespace {
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

} // namespace rhbm_gem::core
