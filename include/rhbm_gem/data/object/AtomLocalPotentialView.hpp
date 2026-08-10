#pragma once

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
    AtomLocalPotentialView() = default;
    static AtomLocalPotentialView For(const AtomObject & atom_object);
    static AtomLocalPotentialView RequireFor(const AtomObject & atom_object);
    bool IsAvailable() const;
    const LocalGaussianResult & GetGaussianResult(LocalFittingStage stage) const;
    const GaussianModel3D & GetEstimateOLS(LocalFittingStage stage) const;
    const GaussianModel3D & GetEstimateMDPDE(LocalFittingStage stage) const;
    LocalPotentialSampleList GetRawSamplingEntries(bool apply_selection = true) const;
    LocalPotentialSampleList GetPeelingSamplingEntries(bool apply_selection = true) const;
    int GetNeighborCountForPeeling() const;
    double GetAlphaR(LocalFittingStage stage) const;

private:
    explicit AtomLocalPotentialView(const AtomObject * atom_object);
    const LocalPotentialEntry * FindEntry() const;
    const LocalPotentialEntry & RequireEntry(const char * context) const;

};

} // namespace rhbm_gem
