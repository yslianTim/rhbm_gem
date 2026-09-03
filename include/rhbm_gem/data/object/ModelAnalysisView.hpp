#pragma once

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
    bool HasGroupedAnalysisData(FittingStage stage) const;
    bool HasAtomGroup(FittingStage stage, GroupKey group_key) const;
    const GaussianModel3D & GetAtomGroupMean(FittingStage stage, GroupKey group_key) const;
    const GaussianModel3D & GetAtomGroupMDPDE(FittingStage stage, GroupKey group_key) const;
    const GaussianModel3D & GetAtomGroupPrior(FittingStage stage, GroupKey group_key) const;
    GaussianModel3DWithUncertainty GetAtomGroupPriorWithUncertainty(
        FittingStage stage, GroupKey group_key) const;
    const std::vector<AtomObject *> & GetAtomObjectList(FittingStage stage, GroupKey group_key) const;
    double GetAtomAlphaR(FittingStage stage, GroupKey group_key) const;
    double GetAtomAlphaG(FittingStage stage, GroupKey group_key) const;
    std::vector<GroupKey> CollectAtomGroupKeys(FittingStage stage) const;

private:
    explicit ModelAnalysisView(const ModelObject & model_object);
    friend class ModelObject;
};

} // namespace rhbm_gem
