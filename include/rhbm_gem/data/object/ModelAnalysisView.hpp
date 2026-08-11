#pragma once

#include <string>
#include <vector>

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem {

class AtomObject;
class ModelObject;

class ModelAnalysisView
{
    const ModelObject & m_model_object;

public:
    explicit ModelAnalysisView(const ModelObject & model_object);
    bool HasGroupedAnalysisData(LocalFittingStage stage) const;
    bool HasAtomGroup(LocalFittingStage stage, GroupKey group_key) const;
    const GaussianModel3D & GetAtomGroupMean(
        LocalFittingStage stage, GroupKey group_key) const;
    const GaussianModel3D & GetAtomGroupMDPDE(
        LocalFittingStage stage, GroupKey group_key) const;
    const GaussianModel3D & GetAtomGroupPrior(
        LocalFittingStage stage, GroupKey group_key) const;
    GaussianModel3DWithUncertainty GetAtomGroupPriorWithUncertainty(
        LocalFittingStage stage, GroupKey group_key) const;
    const std::vector<AtomObject *> & GetAtomObjectList(
        LocalFittingStage stage, GroupKey group_key) const;
    double GetAtomAlphaR(LocalFittingStage stage, GroupKey group_key) const;
    double GetAtomAlphaG(LocalFittingStage stage, GroupKey group_key) const;
    std::vector<GroupKey> CollectAtomGroupKeys(LocalFittingStage stage) const;
    std::string GetAtomCountingSummary() const;
    std::string GetAtomGroupingSummary(LocalFittingStage stage) const;
    
};

} // namespace rhbm_gem
