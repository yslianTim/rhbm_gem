#pragma once

#include <optional>
#include <string>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class LocalPotentialEntry;
class ModelObject;

class AtomLocalPotentialView
{
    const AtomObject * m_atom_object{ nullptr };

public:
    AtomLocalPotentialView() = default;
    static AtomLocalPotentialView For(const AtomObject & atom_object);
    static AtomLocalPotentialView RequireFor(const AtomObject & atom_object);
    bool IsAvailable() const;
    const LocalGaussianResult & GetGaussianResult() const;
    const GaussianModel3D & GetEstimateOLS() const;
    const GaussianModel3D & GetEstimateMDPDE() const;
    const LocalPotentialSampleList & GetSamplingEntries() const;
    double GetAlphaR() const;
    std::optional<LocalPotentialAnnotation> FindAnnotation(const std::string & key) const;

private:
    explicit AtomLocalPotentialView(const AtomObject * atom_object);
    const LocalPotentialEntry * FindEntry() const;
    const LocalPotentialEntry & RequireEntry(const char * context) const;

};

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
