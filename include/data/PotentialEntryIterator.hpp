#pragma once

#include <cstddef>
#include <tuple>
#include <vector>
#include <memory>
#include <map>
#include <string>
#include <unordered_map>

#include "GlobalEnumClass.hpp"
#include "GroupPotentialEntry.hpp"

class AtomObject;
class BondObject;
class ModelObject;
class LocalPotentialEntry;

#ifdef HAVE_ROOT
class TH1D;
class TH2D;
class TGraphErrors;
class TGraph2DErrors;
class TF1;
#endif

class PotentialEntryIterator
{
    AtomObject * m_atom_object;
    BondObject * m_bond_object;
    ModelObject * m_model_object;
    LocalPotentialEntry * m_atom_local_entry;
    LocalPotentialEntry * m_bond_local_entry;

public:
    PotentialEntryIterator(ModelObject * model_object);
    PotentialEntryIterator(AtomObject * atom_object);
    PotentialEntryIterator(BondObject * bond_object);
    ~PotentialEntryIterator();
    double GetAtomGausEstimateMinimum(int par_id, Element element) const;
    double GetBondGausEstimateMinimum(int par_id) const;
    bool IsOutlierAtom(const std::string & class_key) const;
    bool IsOutlierBond(const std::string & class_key) const;
    bool IsAvailableAtomGroupKey(GroupKey group_key, const std::string & class_key, bool varbose=false) const;
    bool IsAvailableBondGroupKey(GroupKey group_key, const std::string & class_key, bool varbose=false) const;
    double GetAtomGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetBondGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetAtomGausVariancePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetBondGausVariancePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    const std::vector<AtomObject *> & GetAtomObjectList(GroupKey group_key, const std::string & class_key) const;
    const std::vector<BondObject *> & GetBondObjectList(GroupKey group_key, const std::string & class_key) const;
    std::vector<AtomObject *> GetOutlierAtomObjectList(GroupKey group_key, const std::string & class_key) const;
    std::unordered_map<int, AtomObject *> GetAtomObjectMap(GroupKey group_key, const std::string & class_key) const;

    std::vector<std::tuple<float, float>> GetLinearModelDistanceAndMapValueList(void) const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList(void) const;
    std::vector<std::tuple<float, float>> GetBinnedDistanceAndMapValueList(int bin_size=15, double x_min=0.0, double x_max=1.5) const;
    std::tuple<float, float> GetDistanceRange(double margin_rate=0.0) const;
    std::tuple<float, float> GetMapValueRange(double margin_rate=0.0) const;
    std::tuple<float, float> GetBinnedMapValueRange(int bin_size=15, double x_min=0.0, double x_max=1.5, double margin_rate=0.0) const;
    double GetAmplitudeEstimateMDPDE(void) const;
    double GetAmplitudeEstimatePosterior(const std::string & key) const;
    double GetAmplitudeVariancePosterior(const std::string & key) const;
    double GetWidthEstimateMDPDE(void) const;
    double GetWidthEstimatePosterior(const std::string & key) const;
    double GetWidthVariancePosterior(const std::string & key) const;
    double GetAlphaR(void) const;
    double GetAtomAlphaG(GroupKey group_key, const std::string & class_key) const;
    double GetBondAlphaG(GroupKey group_key, const std::string & class_key) const;
    Residue GetResidueFromAtomGroupKey(GroupKey group_key, const std::string & class_key) const;
    Residue GetResidueFromBondGroupKey(GroupKey group_key, const std::string & class_key) const;

    #ifdef HAVE_ROOT
    std::unique_ptr<TH1D> CreateComponentCountHistogram(std::vector<GroupKey> & group_key_list, const std::string & class_key) const;
    std::unique_ptr<TH2D> CreateAtomResidueCountHistogram2D(const std::string & class_key);
    std::unique_ptr<TH1D> CreateAtomResidueCountHistogram(const std::string & class_key, Structure structure=static_cast<Structure>(0));
    std::unique_ptr<TH1D> CreateBondResidueCountHistogram(const std::string & class_key);
    std::unique_ptr<TH1D> CreateAtomGausEstimateHistogram(GroupKey group_key, const std::string & class_key, int par_id) const;
    std::unique_ptr<TH1D> CreateLinearModelDataHistogram(int dimension_id) const;
    std::unique_ptr<TH2D> CreateDistanceToMapValueHistogram(int x_bin_size=15, int y_bin_size=1000) const;
    std::vector<std::unique_ptr<TH1D>> CreateMainChainAtomGausRankHistogram(int par_id, int & chain_size, Residue residue=Residue::UNK, size_t extra_id=0, std::vector<Residue> veto_residues_list={});
    std::unique_ptr<TGraphErrors> CreateNormalizedAtomGausEstimateScatterGraph(Element element, double reference_amplitude, bool reverse=false);
    std::unique_ptr<TGraphErrors> CreateNormalizedBondGausEstimateScatterGraph(Element element, double reference_amplitude, bool reverse=false);
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> CreateAtomMapValueToSequenceIDGraphMap(size_t main_chain_element_id, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> CreateAtomGausEstimateToSequenceIDGraphMap(size_t main_chain_element_id, const int par_id=0, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> CreateBondGausEstimateToSequenceIDGraphMap(size_t main_chain_element_id, const int par_id=0, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> CreateAtomGausEstimatePosteriorToSequenceIDGraphMap(size_t main_chain_element_id, const std::string & class_key, const int par_id=0, Residue residue=Residue::UNK);
    std::unique_ptr<TGraphErrors> CreateAtomGausEstimateToResidueGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateBondGausEstimateToResidueGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateAtomGausEstimateToSpotGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateAtomGausEstimateToAtomIdGraph(const std::map<std::string, GroupKey> & group_key_map, const std::vector<std::string> & atom_id_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateAtomGausEstimateScatterGraph(GroupKey group_key, const std::string & class_key, int par1_id=0, int par2_id=1, bool select_outliers=false);
    std::unique_ptr<TGraphErrors> CreateAtomGausEstimateScatterGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, int par1_id=0, int par2_id=1);
    std::unique_ptr<TGraphErrors> CreateBondGausEstimateScatterGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, int par1_id=0, int par2_id=1);
    std::unique_ptr<TGraphErrors> CreateAtomGausEstimateScatterGraph(Element element, bool reverse=false);
    std::unique_ptr<TGraphErrors> CreateBondGausEstimateScatterGraph(Element element, bool reverse=false);
    std::unique_ptr<TGraphErrors> CreateGausEstimateScatterGraph(GroupKey group_key1, GroupKey group_key2, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateDistanceToMapValueGraph(void);
    std::unique_ptr<TGraphErrors> CreateLinearModelDistanceToMapValueGraph(void);
    std::unique_ptr<TGraphErrors> CreateBinnedDistanceToMapValueGraph(int bin_size=15, double x_min=0.0, double x_max=1.5);
    std::unique_ptr<TGraphErrors> CreateInRangeAtomsToGausEstimateGraph(GroupKey group_key, const std::string & class_key, double range=5.0, int par_id=0);
    std::unique_ptr<TGraphErrors> CreateCOMDistanceToGausEstimateGraph(GroupKey group_key, const std::string & class_key, int par_id=0);
    std::unique_ptr<TGraphErrors> CreateAtomXYPositionTomographyGraph(double normalized_z_pos=0.5, double z_ratio_window=0.1, bool com_center=false);
    std::unique_ptr<TF1> CreateAtomLocalLinearModelFunctionOLS(void) const;
    std::unique_ptr<TF1> CreateAtomLocalLinearModelFunctionMDPDE(void) const;
    std::unique_ptr<TF1> CreateAtomLocalGausFunctionOLS(void) const;
    std::unique_ptr<TF1> CreateAtomLocalGausFunctionMDPDE(void) const;
    std::unique_ptr<TF1> CreateAtomGroupLinearModelFunctionMean(GroupKey group_key, const std::string & class_key, double x_min, double x_max) const;
    std::unique_ptr<TF1> CreateAtomGroupLinearModelFunctionPrior(GroupKey group_key, const std::string & class_key, double x_min, double x_max) const;
    std::unique_ptr<TF1> CreateAtomGroupGausFunctionMean(GroupKey group_key, const std::string & class_key) const;
    std::unique_ptr<TF1> CreateAtomGroupGausFunctionPrior(GroupKey group_key, const std::string & class_key) const;
    std::unique_ptr<TF1> CreateBondGroupGausFunctionPrior(GroupKey group_key, const std::string & class_key) const;
    #endif

private:
    size_t GetAtomResidueCount(const std::string & class_key, Residue residue, Structure structure=static_cast<Structure>(0)) const;
    size_t GetBondResidueCount(const std::string & class_key, Residue residue) const;
    bool IsAtomObjectAvailable(void) const;
    bool IsBondObjectAvailable(void) const;
    bool IsAtomLocalEntryAvailable(void) const;
    bool IsBondLocalEntryAvailable(void) const;
    bool IsModelObjectAvailable(void) const;
    bool CheckAtomGroupKey(GroupKey group_key, const std::string & class_key, bool verbose=true) const;
    bool CheckBondGroupKey(GroupKey group_key, const std::string & class_key, bool verbose=true) const;

};
