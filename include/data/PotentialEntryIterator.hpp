#pragma once

#include <cstddef>
#include <tuple>
#include <vector>
#include <memory>
#include <unordered_map>

#include "GlobalEnumClass.hpp"
#include "GroupPotentialEntry.hpp"

class AtomObject;
class BondObject;
class ModelObject;
class LocalPotentialEntry;

#ifdef HAVE_ROOT
class TH1D;
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
    std::unordered_map<int, AtomObject *> GetAtomObjectMap(GroupKey group_key, const std::string & class_key) const;

    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList(void) const;
    std::vector<std::tuple<float, float>> GetBinnedDistanceAndMapValueList(int bin_size=15, double x_min=0.0, double x_max=1.5) const;
    std::tuple<float, float> GetDistanceRange(double margin_rate=0.0) const;
    std::tuple<float, float> GetMapValueRange(double margin_rate=0.0) const;
    double GetAmplitudeEstimateMDPDE(void) const;
    double GetWidthEstimateMDPDE(void) const;

    #ifdef HAVE_ROOT
    std::unique_ptr<TH1D> CreateAtomResidueCountHistogram(const std::string & class_key, Structure structure=static_cast<Structure>(0));
    std::unique_ptr<TH1D> CreateBondResidueCountHistogram(const std::string & class_key);
    std::vector<std::unique_ptr<TH1D>> CreateMainChainRankHistogram(int par_id, int & chain_size, Residue residue=Residue::UNK, size_t extra_id=0, std::vector<Residue> veto_residues_list={});
    std::unique_ptr<TGraphErrors> CreateAmplitudeRatioToWidthScatterGraph(size_t target_id, size_t reference_id, Residue residue=Residue::UNK);
    std::unique_ptr<TGraphErrors> CreateNormalizedGausEstimateScatterGraph(Element element, double reference_amplitude, bool reverse=false);
    std::unique_ptr<TGraphErrors> CreateBfactorToWidthScatterGraph(GroupKey group_key, const std::string & class_key);
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> CreateAtomGausEstimateToResidueIDGraphMap(size_t main_chain_element_id, const int par_id=0, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> CreateBondGausEstimateToResidueIDGraphMap(size_t main_chain_element_id, const int par_id=0, Residue residue=Residue::UNK);
    std::unique_ptr<TGraphErrors> CreateAtomGausEstimateToResidueGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateBondGausEstimateToResidueGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateAtomGausEstimateScatterGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, int par1_id=0, int par2_id=1);
    std::unique_ptr<TGraphErrors> CreateBondGausEstimateScatterGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, int par1_id=0, int par2_id=1);
    std::unique_ptr<TGraphErrors> CreateGausEstimateScatterGraph(Element element, bool reverse=false);
    std::unique_ptr<TGraphErrors> CreateGausEstimateScatterGraph(GroupKey group_key1, GroupKey group_key2, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateDistanceToMapValueGraph(void);
    std::unique_ptr<TGraphErrors> CreateBinnedDistanceToMapValueGraph(int bin_size=15, double x_min=0.0, double x_max=1.5);
    std::unique_ptr<TGraphErrors> CreateInRangeAtomsToGausEstimateGraph(GroupKey group_key, const std::string & class_key, double range=5.0, int par_id=0);
    std::unique_ptr<TGraphErrors> CreateCOMDistanceToGausEstimateGraph(GroupKey group_key, const std::string & class_key, int par_id=0);
    std::unique_ptr<TGraphErrors> CreateXYPositionTomographyGraph(double normalized_z_pos=0.5, double z_ratio_window=0.1, bool com_center=false);
    std::unordered_map<size_t, std::unique_ptr<TGraphErrors>> CreateXYPositionTomographyGraphMap(double normalized_z_pos=0.5, double z_ratio_window=0.1, bool com_center=false);
    std::unique_ptr<TGraphErrors> CreateXYPositionTomographyGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, double normalized_z_pos=0.5, double z_ratio_window=0.1);
    std::unique_ptr<TGraph2DErrors> CreateXYPositionTomographyToGausEstimateGraph2D(std::vector<GroupKey> & group_key_list, const std::string & class_key, double normalized_z_pos=0.5, double z_ratio_window=0.1, int par_id=0);
    std::unordered_map<size_t, std::unique_ptr<TGraph2DErrors>> CreateXYPositionTomographyToGausEstimateGraph2DMap(double normalized_z_pos=0.5, double z_ratio_window=0.1, int par_id=0, bool com_center=false);
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
    Residue GetResidueFromAtomGroupKey(GroupKey group_key, const std::string & class_key) const;
    Residue GetResidueFromBondGroupKey(GroupKey group_key, const std::string & class_key) const;

};
