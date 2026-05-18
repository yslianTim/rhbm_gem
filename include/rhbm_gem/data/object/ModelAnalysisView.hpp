#pragma once

#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

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
    std::tuple<float, float> GetDistanceRange(double margin_rate = 0.0) const;
    std::tuple<float, float> GetResponseRange(double margin_rate = 0.0) const;
    SeriesPointList GetBinnedDistanceResponseSeries(
        int bin_size = 15,
        double x_min = 0.0,
        double x_max = 1.5) const;
    double GetAlphaR() const;
    std::optional<LocalPotentialAnnotation> FindAnnotation(const std::string & key) const;
    double GetMapValueNearCenter() const;
    double CalculateQScore(int par_choice) const;

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
    std::string DescribeAtomGrouping() const;
    double GetAtomGausEstimateMinimum(int par_id, Element element) const;
    bool HasAtomGroup(GroupKey group_key, const std::string & class_key, bool verbose=false) const;
    const GaussianModel3D & GetAtomGroupMean(GroupKey group_key, const std::string & class_key) const;
    const GaussianModel3D & GetAtomGroupMDPDE(GroupKey group_key, const std::string & class_key) const;
    const GaussianModel3D & GetAtomGroupPrior(GroupKey group_key, const std::string & class_key) const;
    GaussianModel3DWithUncertainty GetAtomGroupPriorWithUncertainty(GroupKey group_key, const std::string & class_key) const;
    double GetAtomGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetAtomGausPriorStandardDeviation(GroupKey group_key, const std::string & class_key, int par_id) const;
    const std::vector<AtomObject *> & GetAtomObjectList(GroupKey group_key, const std::string & class_key) const;
    std::vector<AtomObject *> GetOutlierAtomObjectList(GroupKey group_key, const std::string & class_key) const;
    std::unordered_map<int, AtomObject *> GetAtomObjectMap(GroupKey group_key, const std::string & class_key) const;
    double GetAtomAlphaR(GroupKey group_key, const std::string & class_key) const;
    double GetAtomAlphaG(GroupKey group_key, const std::string & class_key) const;
    Residue DecodeResidueFromAtomGroupKey(GroupKey group_key, const std::string & class_key) const;
    Residue GetResidueFromAtomGroupKey(GroupKey group_key, const std::string & class_key) const;
    std::vector<GroupKey> CollectAtomGroupKeys(const std::string & class_key) const;
    
};

} // namespace rhbm_gem
