#pragma once

#include <cstddef>
#include <optional>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem {

class AtomObject;
class LocalPotentialEntry;

class AtomLocalPotentialView
{
    const AtomObject * m_atom_object{ nullptr };

public:
    static AtomLocalPotentialView For(const AtomObject & atom_object);
    bool IsAvailable() const;
    const LocalGaussianResult & GetGaussianResult(FittingStage stage) const;
    const GaussianModel3D & GetEstimateOLS(FittingStage stage) const;
    const GaussianModel3D & GetEstimateMDPDE(FittingStage stage) const;
    LocalPotentialSampleList GetRawSamplingEntries(bool apply_selection = true) const;
    LocalPotentialSampleList GetPeelingSamplingEntries(bool apply_selection = false) const;
    LocalPotentialSampleList GetSamplingEntries(FittingStage stage) const;
    bool HasEnoughSamplingEntriesInRange(
        FittingStage stage,
        double distance_min,
        double distance_max,
        std::size_t minimum_sample_count) const;
    std::optional<double> GetLocalFittingPeelingRatio(
        bool peeling_applied,
        double distance_min,
        double distance_max) const;
    int GetNeighborCountForPeeling() const;
    double GetAlphaR(FittingStage stage) const;

private:
    explicit AtomLocalPotentialView(const AtomObject * atom_object);
    const LocalPotentialEntry * FindEntry() const;
    const LocalPotentialEntry & RequireEntry(const char * context) const;
};

} // namespace rhbm_gem
