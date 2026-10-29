#include <rhbm_gem/core/CommandSystem.hpp>
#include <rhbm_gem/core/CommandTypes.hpp>
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/JointAnalysisFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

template <typename T>
constexpr bool kHeaderExposesCompleteType = sizeof(T) > 0;

int main()
{
    using namespace rhbm_gem;
    using namespace rhbm_gem::core;

    static_assert(kHeaderExposesCompleteType<FitOptions>);
    static_assert(static_cast<int>(PotentialEstimator::TWO_STAGE) == 0);
    static_assert(static_cast<int>(PotentialEstimator::JOINT_COMPONENTS) == 1);

    FitOptions fit_options;
    PotentialAnalysisRequest request;
    if (fit_options.estimator != PotentialEstimator::JOINT_COMPONENTS ||
        request.estimator != PotentialEstimator::JOINT_COMPONENTS)
        return 1;

    using PotentialFittingWorkflow = void (*)(MapObject &, ModelObject &, const FitOptions &);
    PotentialFittingWorkflow volatile workflow_entry{&RunPotentialFittingWorkflow};
    if (workflow_entry == nullptr) return 2;

    JointAnalysisResult result;
    result.atom_ids = {"1"};
    result.row_ids = {"row-0"};
    result.fitting_row_ids = {"row-0"};
    result.fitting_to_original_contributors = {0};
    result.fitting_to_original_rows = {0};
    result.available_row_mask = {true};
    result.selection_domain = JointSelectionDomain{};
    result.selection_domain->target_indices = {0};
    result.selection_domain->contributor_policy =
        std::string(PureHaloSupportExclusionContributorPolicy);
    result.selection_domain->pure_halo_support_exclusion =
        JointPureHaloSupportExclusionProvenance{1, 1, 0, 0, 0, 0};
    result.metadata.software = JointSoftwareProvenance{
        "consumer-smoke", std::string(64, '0'), std::string(64, '0'), std::string(64, '0')};

    WriteJointAnalysisResult(result, "consumer-joint.json", "consumer-joint.csv");
    std::ifstream json("consumer-joint.json");
    const std::string payload((std::istreambuf_iterator<char>(json)), {});
    if (payload.find("\"schema_version\":21") == std::string::npos ||
        payload.find("target-only-pure-halo-support-exclusion-v1") == std::string::npos)
        return 3;

    std::vector<std::unique_ptr<AtomObject>> atoms;
    auto atom = std::make_unique<AtomObject>();
    atom->SetSerialID(1);
    atom->SetElement(Element::CARBON);
    atoms.push_back(std::move(atom));
    ModelObject model(std::move(atoms));
    model.EditAnalysis().ApplyJointResult(result, "consumer-smoke");

    DataRepository repository("consumer-analysis.sqlite");
    repository.SaveModel(model, "consumer");
    const auto loaded = repository.LoadModel("consumer");
    const auto & loaded_result = loaded->GetAnalysisView().GetJointResult();
    if (!loaded_result || loaded_result->atom_ids != result.atom_ids ||
        IsQualifiedJointAnalysisResult(*loaded_result))
        return 4;

    return 0;
}
