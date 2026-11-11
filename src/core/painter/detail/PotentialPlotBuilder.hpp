#pragma once

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
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

enum class DistanceSamplingMode
{
    SelectedRaw,
    AllRaw,
    Peeling
};

enum class QScoreMode { Standard, JointPeeledStandard, JointPeeledModel };
inline constexpr std::array<const char *,3> QScoreMethodLabels{
    "Standard QScore", "Joint-Peeled QScore (standard reference)",
    "Joint-Peeled QScore (Joint model reference)"};
using ResidueQScoreMap = std::map<std::string,std::map<int,double>>;
ResidueQScoreMap CollectResidueQScores(const ModelObject &, QScoreMode);
std::pair<double,double> QScoreDisplayRange(const std::vector<double> &);

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
    std::vector<std::unique_ptr<::TH1D>> CreateMainChainAtomGausRankHistogram(int par_id, int & chain_size, Residue residue=Residue::UNK, size_t extra_id=0, std::vector<Residue> veto_residues_list={});
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAtomMapValueToSequenceIDGraphMap(size_t main_chain_element_id, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAtomGausEstimateToSequenceIDGraphMap(size_t main_chain_element_id, const int par_id=0, Residue residue=Residue::UNK);
    std::unordered_map<std::string, std::unique_ptr<::TGraphErrors>> CreateAverageQScoreToSequenceIDGraphMap(QScoreMode mode);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateToResidueGraph(const std::vector<GroupKey> & group_key_list, const int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateToAtomIdGraph(const std::map<std::string, GroupKey> & group_key_map, const std::vector<std::string> & atom_id_list, const int par_id=0);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateScatterGraph(GroupKey group_key, int par1_id=0, int par2_id=1, bool select_outliers=false);
    std::unique_ptr<::TGraphErrors> CreateAtomGausEstimateScatterGraph(const std::vector<GroupKey> & group_key_list, int par1_id=0, int par2_id=1);
    std::unique_ptr<::TGraphErrors> CreateJointStageParameterScatterGraph(Element element, int choice=0);
    std::unique_ptr<::TGraphErrors> CreateDistanceToMapValueGraph(
        DistanceSamplingMode mode = DistanceSamplingMode::SelectedRaw);
    std::unique_ptr<::TGraphErrors> CreateBinnedDistanceToRawMapValueGraph(int bin_size=20, double x_min=0.0, double x_max=2.0);
    std::unique_ptr<::TGraphErrors> CreateBinnedDistanceToPeelingMapValueGraph(int bin_size=20, double x_min=0.0, double x_max=2.0);
    std::unique_ptr<::TGraphErrors> CreateBinnedDistanceToNeighborMapValueGraph(int bin_size=20, double x_min=0.0, double x_max=2.0);
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
    std::unique_ptr<::TF1> CreateAtomLocalGausFunctionFinal() const;
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
