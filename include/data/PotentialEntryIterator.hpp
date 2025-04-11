#pragma once

#include <tuple>
#include <vector>
#include <memory>
#include <unordered_map>
#include "AtomicInfoHelper.hpp"
#include "GroupPotentialEntry.hpp"

class AtomObject;
class ModelObject;
class AtomicPotentialEntry;

#ifdef HAVE_ROOT
class TGraphErrors;
class TGraph2DErrors;
class TF1;
#endif

class PotentialEntryIterator
{
    using ElementKeyType = GroupKeyMapping<ElementGroupClassifierTag>::type;
    using ResidueKeyType = GroupKeyMapping<ResidueGroupClassifierTag>::type;

    AtomObject * m_atom_object;
    ModelObject * m_model_object;
    AtomicPotentialEntry * m_atomic_entry;
    GroupPotentialEntry<ElementKeyType> * m_element_class_group_entry;
    GroupPotentialEntry<ResidueKeyType> * m_residue_class_group_entry;

public:
    PotentialEntryIterator(ModelObject * model_object);
    PotentialEntryIterator(AtomObject * atom_object);
    ~PotentialEntryIterator();
    bool IsAvailableGroupKey(ElementKeyType & group_key) const;
    bool IsAvailableGroupKey(ResidueKeyType & group_key) const;
    std::tuple<double, double> GetGausEstimatePrior(ElementKeyType & group_key) const;
    std::tuple<double, double> GetGausEstimatePrior(ResidueKeyType & group_key) const;
    std::tuple<double, double> GetGausVariancePrior(ElementKeyType & group_key) const;
    std::tuple<double, double> GetGausVariancePrior(ResidueKeyType & group_key) const;
    double GetGausEstimatePrior(ElementKeyType & group_key, int par_id) const;
    double GetGausEstimatePrior(ResidueKeyType & group_key, int par_id) const;
    double GetGausVariancePrior(ElementKeyType & group_key, int par_id) const;
    double GetGausVariancePrior(ResidueKeyType & group_key, int par_id) const;
    const std::vector<AtomObject *> & GetAtomObjectList(ElementKeyType & group_key) const;
    const std::vector<AtomObject *> & GetAtomObjectList(ResidueKeyType & group_key) const;

    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList(void) const;
    std::tuple<float, float> GetDistanceRange(double margin_rate=0.0) const;
    std::tuple<float, float> GetMapValueRange(double margin_rate=0.0) const;
    double GetAmplitudeEstimateMDPDE(void) const;
    double GetWidthEstimateMDPDE(void) const;

    #ifdef HAVE_ROOT
    std::unique_ptr<TGraphErrors> CreateGausEstimateToResidueGraph(std::vector<ResidueKeyType> & group_key_list, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateGausEstimateScatterGraph(std::vector<ResidueKeyType> & group_key_list, bool reverse=false);
    std::unique_ptr<TGraphErrors> CreateDistanceToMapValueGraph(void);
    std::unique_ptr<TGraphErrors> CreateBinnedDistanceToMapValueGraph(int bin_size=15, double x_min=0.0, double x_max=1.5);
    std::unique_ptr<TGraphErrors> CreateInRangeAtomsToGausEstimateGraph(ElementKeyType & group_key, double range=5.0, int par_id=0);
    std::unique_ptr<TGraphErrors> CreateInRangeAtomsToGausEstimateGraph(ResidueKeyType & group_key, double range=5.0, int par_id=0);
    std::unique_ptr<TGraphErrors> CreateXYPositionTomographyGraph(double normalized_z_pos=0.5, double z_ratio_window=0.1);
    std::unique_ptr<TGraphErrors> CreateXYPositionTomographyGraph(std::vector<ResidueKeyType> & group_key_list, double normalized_z_pos=0.5, double z_ratio_window=0.1);
    std::unique_ptr<TGraph2DErrors> CreateXYPositionTomographyToGausEstimateGraph2D(std::vector<ResidueKeyType> & group_key_list, double normalized_z_pos=0.5, double z_ratio_window=0.1, int par_id=0);
    std::unique_ptr<TF1> CreateGroupGausFunctionPrior(ElementKeyType & group_key) const;
    std::unique_ptr<TF1> CreateGroupGausFunctionPrior(ResidueKeyType & group_key) const;
    #endif

private:
    bool IsAtomObjectAvailable(void) const;
    bool IsAtomicEntryAvailable(void) const;
    bool IsModelObjectAvailable(void) const;
    bool CheckParameterIndex(int par_id) const;
    bool CheckGroupKey(ElementKeyType & group_key, bool verbose=true) const;
    bool CheckGroupKey(ResidueKeyType & group_key, bool verbose=true) const;

};