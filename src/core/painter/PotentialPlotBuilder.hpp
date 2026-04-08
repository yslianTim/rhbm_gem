#pragma once

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#ifdef HAVE_ROOT
class TH1D;
class TH2D;
class TGraphErrors;
class TF1;
#endif

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ModelObject;

class PotentialPlotBuilder
{
    ModelObject * m_model_object{ nullptr };
    AtomObject * m_atom_object{ nullptr };
    BondObject * m_bond_object{ nullptr };

public:
    explicit PotentialPlotBuilder(ModelObject * model_object);
    explicit PotentialPlotBuilder(AtomObject * atom_object);
    explicit PotentialPlotBuilder(BondObject * bond_object);

#ifdef HAVE_ROOT
    std::unique_ptr<::TH1D> CreateComponentCountHistogram(std::vector<GroupKey> & group_key_list, const std::string & class_key) const;
    std::unique_ptr<::TH2D> CreateAtomResidueCountHistogram2D(const std::string & class_key);
    std::unique_ptr<::TH1D> CreateAtomResidueCountHistogram(const std::string & class_key, Structure structure=static_cast<Structure>(0));
    std::unique_ptr<::TH1D> CreateBondResidueCountHistogram(const std::string & class_key);
    std::unique_ptr<::TH1D> CreateAtomGausEstimateHistogram(GroupKey group_key, const std::string & class_key, int par_id) const;
    std::unique_ptr<::TH1D> CreateLinearModelDataHistogram(int dimension_id) const;
    std::unique_ptr<::TH2D> CreateDistanceToMapValueHistogram(int x_bin_size=15, int y_bin_size=1000) const;
    std::vector<std::unique_ptr<::TH1D>> CreateMainChainAtomGausRankHistogram(int par_id, int & chain_size, Residue residue=Residue::UNK, size_t extra_id=0, std::vector<Residue> veto_residues_list={});
    std::unique_ptr<::TGraphErrors> CreateNormalizedAtomGausEstimateScatterGraph(Element element, double reference_amplitude, bool reverse=false);
    std::unique_ptr<::TGraphErrors> CreateNormalizedBondGausEstimateScatterGraph(Element element, double reference_amplitude, bool reverse=false);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAtomMapValueToSequenceIDGraphMap(size_t main_chain_element_id, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAtomGausEstimateToSequenceIDGraphMap(size_t main_chain_element_id, const int par_id=0, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAtomQScoreToSequenceIDGraphMap(size_t main_chain_element_id, const int par_choice=0);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateBondGausEstimateToSequenceIDGraphMap(size_t main_chain_element_id, const int par_id=0, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAtomGausEstimatePosteriorToSequenceIDGraphMap(size_t main_chain_element_id, const std::string & class_key, const int par_id=0, Residue residue=Residue::UNK);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateToResidueGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateBondGausEstimateToResidueGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateToSpotGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateToAtomIdGraph(const std::map<std::string, GroupKey> & group_key_map, const std::vector<std::string> & atom_id_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateScatterGraph(GroupKey group_key, const std::string & class_key, int par1_id=0, int par2_id=1, bool select_outliers=false);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateScatterGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, int par1_id=0, int par2_id=1);
    std::unique_ptr<::TGraphErrors> CreateBondGausEstimateScatterGraph(std::vector<GroupKey> & group_key_list, const std::string & class_key, int par1_id=0, int par2_id=1);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateScatterGraph(Element element, bool reverse=false);
    std::unique_ptr<::TGraphErrors> CreateBondGausEstimateScatterGraph(Element element, bool reverse=false);
    std::unique_ptr<::TGraphErrors> CreateGausEstimateScatterGraph(GroupKey group_key1, GroupKey group_key2, const std::string & class_key, const int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateDistanceToMapValueGraph();
    std::unique_ptr<::TGraphErrors> CreateLinearModelDistanceToMapValueGraph();
    std::unique_ptr<::TGraphErrors> CreateBinnedDistanceToMapValueGraph(int bin_size=15, double x_min=0.0, double x_max=1.5);
    std::unique_ptr<::TGraphErrors> CreateInRangeAtomsToGausEstimateGraph(GroupKey group_key, const std::string & class_key, double range=5.0, int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateCOMDistanceToGausEstimateGraph(GroupKey group_key, const std::string & class_key, int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateAtomXYPositionTomographyGraph(double normalized_z_pos=0.5, double z_ratio_window=0.1, bool com_center=false);
    std::unique_ptr<::TF1> CreateAtomLocalLinearModelFunctionOLS() const;
    std::unique_ptr<::TF1> CreateAtomLocalLinearModelFunctionMDPDE() const;
    std::unique_ptr<::TF1> CreateAtomLocalGausFunctionOLS() const;
    std::unique_ptr<::TF1> CreateAtomLocalGausFunctionMDPDE() const;
    std::unique_ptr<::TF1> CreateAtomGroupLinearModelFunctionMean(GroupKey group_key, const std::string & class_key, double x_min, double x_max) const;
    std::unique_ptr<::TF1> CreateAtomGroupLinearModelFunctionPrior(GroupKey group_key, const std::string & class_key, double x_min, double x_max) const;
    std::unique_ptr<::TF1> CreateAtomGroupGausFunctionMean(GroupKey group_key, const std::string & class_key) const;
    std::unique_ptr<::TF1> CreateAtomGroupGausFunctionPrior(GroupKey group_key, const std::string & class_key) const;
    std::unique_ptr<::TF1> CreateBondGroupGausFunctionPrior(GroupKey group_key, const std::string & class_key) const;
#endif

private:
    ModelAnalysisView GetModelView() const;
    LocalPotentialView GetLocalEntry() const;
    bool IsModelObjectAvailable() const;
    bool IsAtomLocalEntryAvailable() const;
    size_t GetAtomResidueCount(const std::string & class_key, Residue residue, Structure structure=static_cast<Structure>(0)) const;
    size_t GetBondResidueCount(const std::string & class_key, Residue residue) const;
    bool IsAvailableAtomGroupKey(GroupKey group_key, const std::string & class_key, bool varbose=false) const;
    bool IsAvailableBondGroupKey(GroupKey group_key, const std::string & class_key, bool varbose=false) const;
};

} // namespace rhbm_gem
