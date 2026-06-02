#include "detail/CommandBase.hpp"

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

#include <memory>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>

namespace rhbm_gem::core {

class PotentialAnalysisCommand final : public CommandBase<PotentialAnalysisRequest>
{
public:
    PotentialAnalysisCommand();

private:
    void NormalizeAndValidateRequest(PotentialAnalysisRequest & request) override;
    void ValidatePreparedRequest(const PotentialAnalysisRequest & request) override;
    bool ExecuteImpl(const PotentialAnalysisRequest & request) override;
};

namespace {

struct PotentialAnalysisInputs
{
    std::unique_ptr<ModelObject> model_object;
    std::unique_ptr<MapObject> map_object;
};

std::optional<PotentialAnalysisInputs> LoadPotentialAnalysisInputs(
    const PotentialAnalysisRequest & request)
{
    try
    {
        auto model_object{ ReadModel(request.model_file_path) };
        model_object->SetKeyTag("model");
        auto map_object{ ReadMap(request.map_file_path) };
        map_object->SetKeyTag("map");
        if (model_object == nullptr || map_object == nullptr)
        {
            Logger::Log(LogLevel::Error,
                "LoadPotentialAnalysisInputs : model/map object missing after load.");
            return std::nullopt;
        }
        if (request.simulation_flag)
        {
            if (request.simulated_map_resolution == 0.0)
            {
                Logger::Log(LogLevel::Warning,
                    "[Warning] The resolution of input simulated map hasn't been set.\n"
                    "          Please give the corresponding resolution value for this map.\n"
                    "          (-r, --sim-resolution)");
            }
            model_object->SetEmdID("Simulation");
            model_object->SetResolution(request.simulated_map_resolution);
            model_object->SetResolutionMethod("Blurring Width");
        }
        return PotentialAnalysisInputs{ std::move(model_object), std::move(map_object) };
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error, "LoadPotentialAnalysisInputs : " + std::string(e.what()));
        return std::nullopt;
    }
}

void RunModelObjectPreprocessing(
    ModelObject & model_object,
    bool asymmetry_flag,
    bool exclude_hydrogen,
    double alpha_r,
    double alpha_g)
{
    auto analysis{ model_object.EditAnalysis() };
    analysis.Clear();
    model_object.SelectAllAtoms();
    model_object.ApplySymmetrySelection(asymmetry_flag);
    if (exclude_hydrogen)
    {
        const auto selected_atoms{ model_object.GetSelectedAtoms() };
        std::unordered_set<const AtomObject *> selected_atom_set{
            selected_atoms.begin(),
            selected_atoms.end()
        };
        model_object.SelectAtoms(
            [&selected_atom_set](const AtomObject & atom)
            {
                return selected_atom_set.find(&atom) != selected_atom_set.end()
                    && atom.GetElement() != Element::HYDROGEN;
            });
    }

    for (auto * atom : model_object.GetSelectedAtoms())
    {
        analysis.EnsureAtomLocalPotential(*atom);
    }

    // Establish the model-analysis preprocessing invariant for downstream steps:
    // selection is finalized, local entries exist, atom groups are materialized,
    // and selected atoms carry the initial alpha-r and alpha_g.
    analysis.RebuildAtomGroupsFromSelection();
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        analysis.EnsureAtomLocalPotential(*atom).SetAlphaR(alpha_r);
    }
    const auto analysis_view{ model_object.GetAnalysisView() };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        for (const auto group_key : analysis_view.CollectAtomGroupKeys(class_key))
        {
            analysis.SetAtomGroupAlphaG(group_key, class_key, alpha_g);
        }
    }

    Logger::Log(LogLevel::Info,
        "Number of selected atom = " + std::to_string(model_object.GetSelectedAtomCount()));
    std::string description{ "Atom Grouping Summary:" };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        description +=
            "\n - Class type: " + class_key + " include "
            + std::to_string(analysis_view.CollectAtomGroupKeys(class_key).size()) + " groups.";
    }
    Logger::Log(LogLevel::Info, description);
    if (model_object.GetNumberOfAtom() > 0 && model_object.GetSelectedAtomCount() == 0)
    {
        Logger::Log(LogLevel::Warning,
            "No atoms are selected after symmetry filtering. "
            "The input CIF may miss usable _entity/_struct_asym metadata. "
            "Try '--asymmetry true' to bypass symmetry filtering.");
    }
}

} // namespace

PotentialAnalysisCommand::PotentialAnalysisCommand() : CommandBase<PotentialAnalysisRequest>{}
{
}

void PotentialAnalysisCommand::NormalizeAndValidateRequest(PotentialAnalysisRequest & request)
{
    RequireExistingPath(request, &PotentialAnalysisRequest::model_file_path);
    RequireExistingPath(request, &PotentialAnalysisRequest::map_file_path);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::simulated_map_resolution);
    RequireNonEmptyList(request, &PotentialAnalysisRequest::saved_key_tag);
    RequireEnum(request, &PotentialAnalysisRequest::sampling_method);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::fit_range_min);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::fit_range_max);
    RequireFinitePositiveScalar(request, &PotentialAnalysisRequest::alpha_r);
    RequireFinitePositiveScalar(request, &PotentialAnalysisRequest::alpha_g);
}

bool PotentialAnalysisCommand::ExecuteImpl(const PotentialAnalysisRequest & request)
{
    auto inputs{ LoadPotentialAnalysisInputs(request) };
    if (!inputs.has_value()) return false;

    auto & model_object{ *inputs->model_object };
    auto & map_object{ *inputs->map_object };
    if (!request.simulation_flag && request.map_normalization_flag)
    {
        map_object.MapValueArrayNormalization();
    }
    RunModelObjectPreprocessing(
        model_object,
        request.asymmetry_flag,
        request.exclude_hydrogen,
        request.alpha_r,
        request.alpha_g);
    RunPotentialSamplingWorkflow(map_object, model_object, request.sampling_method, request.job_count);

    FitOptions options;
    options.local_fit_model = LocalGaussianFitModel::LogQuadratic;
    options.distance_min = request.fit_range_min;
    options.distance_max = request.fit_range_max;
    options.thread_size = request.job_count;
    if (request.training_alpha_flag)
    {
        RunLocalAlphaTraining(model_object, options);
        RunGroupAlphaTraining(model_object, options);
    }
    else
    {
        RunLocalPotentialFitting(model_object, options);
    }
    RunGroupPotentialFitting(model_object, options);

    DataRepository repository{ request.database_path };
    repository.SaveModel(model_object, request.saved_key_tag);
    model_object.EditAnalysis().ClearTransientFitStates();
    return true;
}

void PotentialAnalysisCommand::ValidatePreparedRequest(const PotentialAnalysisRequest & request)
{
    RequirePrepareCondition(
        !request.simulation_flag || request.simulated_map_resolution > 0.0,
        "Expected a positive simulated-map resolution when '--simulation true' is selected.");
    RequirePrepareCondition(
        request.fit_range_min <= request.fit_range_max,
        "Expected --fit-min <= --fit-max.");
}

namespace command_internal {

CommandResult ExecutePotentialAnalysisCommand(const PotentialAnalysisRequest & request)
{
    PotentialAnalysisCommand command;
    return command.ExecuteRequest(request);
}

} // namespace command_internal

} // namespace rhbm_gem::core
