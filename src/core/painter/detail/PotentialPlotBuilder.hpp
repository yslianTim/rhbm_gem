#pragma once

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#ifdef HAVE_ROOT
class TH1D;
class TH2D;
class TGraphErrors;
class TF1;
#endif

namespace rhbm_gem {

class AtomObject;
class ModelObject;

class PotentialPlotBuilder
{
    ModelObject * m_model_object{ nullptr };
    AtomObject * m_atom_object{ nullptr };

public:
    explicit PotentialPlotBuilder(ModelObject * model_object);
    explicit PotentialPlotBuilder(AtomObject * atom_object);

#ifdef HAVE_ROOT
    std::unique_ptr<::TH1D> CreateComponentCountHistogram(const std::vector<GroupKey> & group_key_list) const;
    std::unique_ptr<::TH1D> CreateAtomGausEstimateHistogram(GroupKey group_key, int par_id) const;
    std::unique_ptr<::TH1D> CreateLinearModelDataHistogram(int dimension_id, bool apply_selection=true, bool use_updated_sample=false) const;
    std::unique_ptr<::TH2D> CreateDistanceToMapValueHistogram(int x_bin_size=20, int y_bin_size=1000, bool apply_selection=true, bool use_updated_sample=false) const;
    std::vector<std::unique_ptr<::TH1D>> CreateMainChainAtomGausRankHistogram(int par_id, int & chain_size, Residue residue=Residue::UNK, size_t extra_id=0, std::vector<Residue> veto_residues_list={});
    std::unique_ptr<::TGraphErrors> CreateNormalizedAtomGausEstimateScatterGraph(Element element, double reference_amplitude, bool reverse=false);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAtomMapValueToSequenceIDGraphMap(size_t main_chain_element_id, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAtomGausEstimateToSequenceIDGraphMap(size_t main_chain_element_id, const int par_id=0, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAtomQScoreToSequenceIDGraphMap(size_t main_chain_element_id, const int par_choice=0, bool apply_selection=true, bool use_updated_sample=false);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAverageQScoreToSequenceIDGraphMap(bool use_fitted_par=false, bool apply_selection=true, bool use_updated_sample=false);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateToResidueGraph(const std::vector<GroupKey> & group_key_list, const int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateToAtomIdGraph(const std::map<std::string, GroupKey> & group_key_map, const std::vector<std::string> & atom_id_list, const int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateScatterGraph(GroupKey group_key, int par1_id=0, int par2_id=1, bool select_outliers=false);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateScatterGraph(const std::vector<GroupKey> & group_key_list, int par1_id=0, int par2_id=1);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateScatterGraph(Element element, bool reverse=false);
    std::unique_ptr<::TGraphErrors> CreateDistanceToMapValueGraph(bool apply_selection=true, bool use_updated_sample=false);
    std::unique_ptr<::TGraphErrors> CreateLinearModelDistanceToMapValueGraph(bool apply_selection=true, bool use_updated_sample=false);
    std::unique_ptr<::TGraphErrors> CreateBinnedDistanceToMapValueGraph(int bin_size=20, double x_min=0.0, double x_max=2.0);
    std::unique_ptr<::TGraphErrors> CreateInRangeAtomsToGausEstimateGraph(GroupKey group_key, double range=5.0, int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateCOMDistanceToGausEstimateGraph(GroupKey group_key, int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateAtomXYPositionTomographyGraph(double normalized_z_pos=0.5, double z_ratio_window=0.1, bool com_center=false);
    static std::unique_ptr<::TGraphErrors> CreateMapValueScatterGraph(
        AtomKey atom_key,
        ModelObject * model1,
        ModelObject * model2,
        int bin_size=15,
        double x_min=0.0,
        double x_max=1.5);
    static std::vector<AtomObject *> CollectComponentAtomMembers(
        const ModelAnalysisView & model_view,
        AtomKey atom_key);
    static std::optional<GaussianModel3DWithUncertainty> ComputeComponentAtomAveragePrior(
        const ModelAnalysisView & model_view,
        AtomKey atom_key);
    std::unique_ptr<::TF1> CreateAtomLocalLinearModelFunctionOLS() const;
    std::unique_ptr<::TF1> CreateAtomLocalLinearModelFunctionMDPDE() const;
    std::unique_ptr<::TF1> CreateAtomLocalGausFunctionOLS() const;
    std::unique_ptr<::TF1> CreateAtomLocalGausFunctionMDPDE() const;
    std::unique_ptr<::TF1> CreateAtomGroupGausFunctionMean(GroupKey group_key) const;
    std::unique_ptr<::TF1> CreateAtomGroupGausFunctionPrior(GroupKey group_key) const;
    std::unique_ptr<::TF1> CreateComponentAtomAverageGausFunctionPrior(AtomKey atom_key) const;
#endif

private:
    ModelAnalysisView GetModelView() const;
    AtomLocalPotentialView GetLocalEntry() const;
    bool IsModelObjectAvailable() const;
    bool IsAtomLocalEntryAvailable() const;
    bool IsAvailableAtomGroupKey(GroupKey group_key) const;
};

} // namespace rhbm_gem
