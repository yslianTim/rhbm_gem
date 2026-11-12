#include <rhbm_gem/core/CommandSystem.hpp>
#include <rhbm_gem/core/CommandTypes.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/JointAnalysisFileIO.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
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

int main()
{
    using namespace rhbm_gem;
    using namespace rhbm_gem::core;

    PotentialAnalysisRequest request;
    request.verbosity = 0;
    const auto rejected = RunCommand(request);
    if (rejected.succeeded || rejected.issues.empty()) return 2;

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
    result.metadata.solver = JointSolverProvenance{};
    auto & solver = *result.metadata.solver;
    solver.search_method = "gauss-seidel-bcd";
    solver.contract_version = std::string(JointSolverProvenanceContractVersion);
    solver.sparse_backend = "EIGEN";
    solver.bcd_core_atoms = 12;
    solver.bcd_policy_version = std::string(BcdPolicyContractVersion);
    solver.bcd_maximum_sweeps = 300;
    solver.bcd_order = "forward";
    solver.bcd_local_search = std::string(BcdBoundedLocalSearchName);

    WriteJointAnalysisResult(result, "consumer-joint.json", "consumer-joint.csv");
    std::ifstream json("consumer-joint.json");
    const std::string payload((std::istreambuf_iterator<char>(json)), {});
    if (payload.find("\"schema_version\":25") == std::string::npos ||
        payload.find("\"search_method\":\"gauss-seidel-bcd\"") == std::string::npos ||
        payload.find("joint-solver-provenance-v6") == std::string::npos ||
        payload.find("target-only-pure-halo-support-exclusion-v1") == std::string::npos)
        return 3;

    std::vector<std::unique_ptr<AtomObject>> atoms;
    auto atom = std::make_unique<AtomObject>();
    atom->SetSerialID(1);
    atom->SetElement(Element::CARBON);
    atoms.push_back(std::move(atom));
    ModelObject model(std::move(atoms));
    model.EditAnalysis().ApplyJointResult(result, "consumer-smoke");
    const auto local = AtomLocalPotentialView::For(*model.FindAtomPtr(1));
    model.EditAnalysis().SetAtomLocalSamplingEntries(*model.FindAtomPtr(1), {
        {6.0, {0.1, {1.0, 2.0, 3.0}, true}},
        {4.0, {0.2, {4.0, 5.0, 6.0}, false}}});
    if (local.GetSelectedSamplingEntries().size() != 1 ||
        local.GetSamplingEntriesView().size() != 2 ||
        !local.GetPeeledObservationEntries().empty())
        return 6;
    const auto & source = local.GetJointStageEstimate().source;
    if (source.method != EstimateMethod::JointComponents || source.atom_id != "1" ||
        source.run_id != "consumer-smoke" || source.role != FittingRole::Target ||
        local.GetInitialStageEstimate().source.method != EstimateMethod::Unspecified)
        return 5;

    DataRepository repository("consumer-analysis.sqlite");
    repository.SaveModel(model, "consumer");
    const auto loaded = repository.LoadModel("consumer");
    const auto & loaded_result = loaded->GetAnalysisView().GetJointResult();
    const auto & loaded_source =
        AtomLocalPotentialView::For(*loaded->FindAtomPtr(1)).GetJointStageEstimate().source;
    if (!loaded_result || loaded_result->atom_ids != result.atom_ids ||
        !loaded_result->metadata.solver ||
        loaded_result->metadata.solver->contract_version != JointSolverProvenanceContractVersion ||
        loaded_result->metadata.solver->bcd_policy_version != BcdPolicyContractVersion ||
        loaded_source.method != source.method || loaded_source.atom_id != source.atom_id ||
        loaded_source.run_id != source.run_id || loaded_source.role != source.role ||
        IsQualifiedJointAnalysisResult(*loaded_result))
        return 4;

    const auto & samples = AtomLocalPotentialView::For(*loaded->FindAtomPtr(1)).GetSamplingEntriesView();
    if (samples.size() != 2 || samples[1].response != 4.0 || samples[1].point.is_selected ||
        samples[1].point.position != local.GetSamplingEntriesView()[1].point.position)
        return 7;

    return 0;
}
