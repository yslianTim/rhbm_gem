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
    bool HasGroupedAnalysisData() const;
    bool HasAtomGroup(GroupKey group_key, const std::string & class_key) const;
    const GaussianModel3D & GetAtomGroupMean(GroupKey group_key, const std::string & class_key) const;
    const GaussianModel3D & GetAtomGroupMDPDE(GroupKey group_key, const std::string & class_key) const;
    const GaussianModel3D & GetAtomGroupPrior(GroupKey group_key, const std::string & class_key) const;
    GaussianModel3DWithUncertainty GetAtomGroupPriorWithUncertainty(GroupKey group_key, const std::string & class_key) const;
    const std::vector<AtomObject *> & GetAtomObjectList(GroupKey group_key, const std::string & class_key) const;
    double GetAtomAlphaR(GroupKey group_key, const std::string & class_key) const;
    double GetAtomAlphaG(GroupKey group_key, const std::string & class_key) const;
    std::vector<GroupKey> CollectAtomGroupKeys(const std::string & class_key) const;
    
};

} // namespace rhbm_gem
