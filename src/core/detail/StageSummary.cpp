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

std::string BuildJointTargetSpotSummary(const ModelObject & model_object)
{
    constexpr int kSpotColumnWidth = 4;
    constexpr int kValidColumnWidth = 5;
    constexpr int kNotConvergedColumnWidth = 13;
    constexpr int kUnavailableColumnWidth = 11;
    constexpr int kMeanColumnWidth = 6;
    constexpr int kStandardDeviationColumnWidth = 4;
    constexpr int kParameterColumnWidth = kMeanColumnWidth + 3 + kStandardDeviationColumnWidth;
    std::map<Spot, GaussianModelParameterSamples> spots;
    for (const auto * atom : model_object.GetSelectedAtoms())
    {
        if (atom->GetSpecialAtomFlag()) continue;
        const auto spot = atom->GetSpot();
        if (std::find(kLocalMDPDESummarySpotList.begin(), kLocalMDPDESummarySpotList.end(), spot) ==
            kLocalMDPDESummarySpotList.end()) continue;
        const auto view = AtomLocalPotentialView::For(*atom);
        const auto * estimate = view.IsAvailable() ? &view.GetStageEstimate(FittingStage::Second) : nullptr;
        const bool joint = estimate && estimate->source.method == EstimateMethod::JointComponents;
        if (joint && estimate->source.role != FittingRole::Target) continue;
        auto & samples = spots[spot];
        if (!joint) { ++samples.unavailable; continue; }
        if (estimate->convergence != JointCheckStatus::Passed || !estimate->point)
        {
            if (estimate->convergence == JointCheckStatus::Failed) ++samples.not_converged;
            else ++samples.unavailable;
            continue;
        }
        const auto & point = *estimate->point;
        samples.amplitude_list.push_back(point.GetAmplitude());
        samples.width_list.push_back(point.GetWidth());
        samples.offset_list.push_back(point.GetOffset());
    }
    std::ostringstream summary;
    summary << "Joint Target estimate summary by Spot:\nEstimator: joint-components"
        "\nPopulation: selected targets excluding special atoms; A/tau/C: Joint Target estimates; "
        "C(J): Joint charge coefficient; s.d.: between-atom dispersion"
        << "\n| " << std::left << std::setw(kSpotColumnWidth) << "Spot"
        << " | " << std::setw(kValidColumnWidth) << "valid"
        << " | " << std::setw(kNotConvergedColumnWidth) << "not-converged"
        << " | " << std::setw(kUnavailableColumnWidth) << "unavailable"
        << " | " << std::setw(kParameterColumnWidth) << "A mean / s.d."
        << " | " << std::setw(kParameterColumnWidth) << "tau mean / s.d."
        << " | " << std::setw(kParameterColumnWidth) << "C(J) mean / s.d."
        << " |";
    for (const auto & [spot, samples] : spots)
    {
        summary << "\n| " << std::left << std::setw(kSpotColumnWidth) << ChemicalDataHelper::GetLabel(spot)
            << " | " << std::right << std::setw(kValidColumnWidth) << samples.amplitude_list.size()
            << " | " << std::setw(kNotConvergedColumnWidth) << samples.not_converged
            << " | " << std::setw(kUnavailableColumnWidth) << samples.unavailable;
        for (const auto * values : {&samples.amplitude_list, &samples.width_list, &samples.offset_list})
        {
            if (values->empty())
            {
                summary << " | " << std::left << std::setw(kParameterColumnWidth) << "unavailable";
                continue;
            }
            const auto mean = array_helper::ComputeMean(values->data(), values->size());
            summary << " | " << std::right << std::fixed << std::setprecision(2)
                << std::setw(kMeanColumnWidth) << mean << " / "
                << std::setw(kStandardDeviationColumnWidth)
                << array_helper::ComputeStandardDeviation(values->data(), values->size(), mean);
        }
        summary << " |";
    }
    if (spots.empty()) summary << "\nNo matching selected targets available.";
    return summary.str();
}

} // namespace rhbm_gem::core
