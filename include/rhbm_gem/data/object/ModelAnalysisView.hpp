#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/data/object/GaussianStatistics.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ModelObject;

struct LocalPotentialAnnotationView
{
    GaussianPosterior posterior{};
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
};

struct GroupingSummaryItem
{
    std::string class_key{};
    size_t group_count{ 0 };
};

struct GroupingSummary
{
    std::string title{};
    std::vector<GroupingSummaryItem> items{};
};

class LocalPotentialView
{
public:
    LocalPotentialView() = default;

    static LocalPotentialView For(const AtomObject & atom_object);
    static LocalPotentialView For(const BondObject & bond_object);
    static LocalPotentialView RequireFor(const AtomObject & atom_object);
    static LocalPotentialView RequireFor(const BondObject & bond_object);

    bool IsAvailable() const;
    const GaussianEstimate & GetEstimateOLS() const;
    const GaussianEstimate & GetEstimateMDPDE() const;
    const LocalPotentialSampleList & GetSamplingEntries() const;
    std::tuple<float, float> GetDistanceRange(double margin_rate = 0.0) const;
    std::tuple<float, float> GetResponseRange(double margin_rate = 0.0) const;
    SeriesPointList GetBinnedDistanceResponseSeries(
        int bin_size = 15,
        double x_min = 0.0,
        double x_max = 1.5) const;
    SeriesPointList GetLinearModelSeries() const;
    int GetSamplingEntryCount() const;
    double GetAlphaR() const;
    std::optional<LocalPotentialAnnotationView> FindAnnotation(const std::string & key) const;
    double GetMapValueNearCenter() const;
    double GetMomentZeroEstimate() const;
    double GetMomentTwoEstimate() const;
    double CalculateQScore(int par_choice) const;
    const AtomObject * GetAtomObjectPtr() const { return m_atom_object; }
    const BondObject * GetBondObjectPtr() const { return m_bond_object; }

private:
    const AtomObject * m_atom_object{ nullptr };
    const BondObject * m_bond_object{ nullptr };

    explicit LocalPotentialView(const AtomObject * atom_object);
    explicit LocalPotentialView(const BondObject * bond_object);
};

class ModelAnalysisView
{
    const ModelObject & m_model_object;

public:
    explicit ModelAnalysisView(const ModelObject & model_object);

    static ModelAnalysisView Of(const ModelObject & model_object);

    bool HasGroupedAnalysisData() const;
    GroupingSummary CollectAtomGroupingSummary() const;
    GroupingSummary CollectBondGroupingSummary() const;
    std::string DescribeAtomGrouping() const;
    std::string DescribeBondGrouping() const;
    double GetAtomGausEstimateMinimum(int par_id, Element element) const;
    double GetBondGausEstimateMinimum(int par_id) const;
    bool HasAtomGroup(GroupKey group_key, const std::string & class_key, bool verbose=false) const;
    bool HasBondGroup(GroupKey group_key, const std::string & class_key, bool verbose=false) const;
    const GaussianEstimate & GetAtomGroupMean(GroupKey group_key, const std::string & class_key) const;
    const GaussianEstimate & GetBondGroupMean(GroupKey group_key, const std::string & class_key) const;
    const GaussianEstimate & GetAtomGroupMDPDE(GroupKey group_key, const std::string & class_key) const;
    const GaussianEstimate & GetBondGroupMDPDE(GroupKey group_key, const std::string & class_key) const;
    const GaussianEstimate & GetAtomGroupPrior(GroupKey group_key, const std::string & class_key) const;
    const GaussianEstimate & GetBondGroupPrior(GroupKey group_key, const std::string & class_key) const;
    GaussianPosterior GetAtomGroupPriorPosterior(GroupKey group_key, const std::string & class_key) const;
    GaussianPosterior GetBondGroupPriorPosterior(GroupKey group_key, const std::string & class_key) const;
    double GetAtomGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetBondGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetAtomGausVariancePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetBondGausVariancePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    const std::vector<AtomObject *> & GetAtomObjectList(GroupKey group_key, const std::string & class_key) const;
    const std::vector<BondObject *> & GetBondObjectList(GroupKey group_key, const std::string & class_key) const;
    std::vector<AtomObject *> GetOutlierAtomObjectList(GroupKey group_key, const std::string & class_key) const;
    std::unordered_map<int, AtomObject *> GetAtomObjectMap(GroupKey group_key, const std::string & class_key) const;
    double GetAtomAlphaR(GroupKey group_key, const std::string & class_key) const;
    double GetAtomAlphaG(GroupKey group_key, const std::string & class_key) const;
    double GetBondAlphaG(GroupKey group_key, const std::string & class_key) const;
    Residue DecodeResidueFromAtomGroupKey(GroupKey group_key, const std::string & class_key) const;
    Residue DecodeResidueFromBondGroupKey(GroupKey group_key, const std::string & class_key) const;
    Residue GetResidueFromAtomGroupKey(GroupKey group_key, const std::string & class_key) const;
    Residue GetResidueFromBondGroupKey(GroupKey group_key, const std::string & class_key) const;
    std::vector<GroupKey> CollectAtomGroupKeys(const std::string & class_key) const;
    std::vector<GroupKey> CollectBondGroupKeys(const std::string & class_key) const;
};

} // namespace rhbm_gem
