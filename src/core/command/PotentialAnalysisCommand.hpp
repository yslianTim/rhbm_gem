#pragma once

#include <memory>

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class ModelObject;
class MapObject;
class AtomObject;
class MutableLocalPotentialView;

class PotentialAnalysisCommand : public CommandWithRequest<PotentialAnalysisRequest>
{
    std::string m_model_key_tag, m_map_key_tag;
    std::shared_ptr<MapObject> m_map_object;
    std::shared_ptr<ModelObject> m_model_object;

public:
    PotentialAnalysisCommand();
    ~PotentialAnalysisCommand() override = default;

private:
    void NormalizeRequest() override;
    void ValidateOptions() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObject(const PotentialAnalysisRequest & request);
    std::vector<MutableLocalPotentialView> BuildSelectedAtomLocalEntryViews(ModelObject & model_object);
    void RunFittingWorkflow(
        ModelObject & model_object,
        MapObject & map_object,
        const PotentialAnalysisRequest & request);
    void UpdateModelObjectForSimulation(
        ModelObject & model_object,
        double simulated_map_resolution);
    void RunMapObjectPreprocessing(MapObject & map_object);
    void RunModelObjectPreprocessing(ModelObject & model_object, bool asymmetry_flag);
    void RunSamplingWorkflow(
        MapObject & map_object,
        ModelObject & model_object,
        int sampling_size,
        double sampling_range_min,
        double sampling_range_max);
    void RunDatasetPreparationWorkflow(
        ModelObject & model_object,
        double fit_range_min,
        double fit_range_max);
    void RunAtomAlphaTraining(
        ModelObject & model_object,
        const std::filesystem::path & training_report_dir);
    double TrainUniversalAlphaR(
        ModelObject & model_object,
        const std::vector<AtomObject *> & atom_list,
        const size_t subset_size,
        const std::vector<double> & alpha_list
    );
    double TrainUniversalAlphaG(
        ModelObject & model_object,
        const std::vector<std::vector<AtomObject *>> & atom_list_set,
        const size_t subset_size,
        const std::vector<double> & alpha_list
    );
    void RunLocalFitting(ModelObject & model_object, double universal_alpha_r);
    void RunAtomPotentialFitting(
        ModelObject & model_object,
        bool training_alpha_flag,
        double fallback_alpha_g);
    void RunExperimentalBondWorkflowIfEnabled(
        ModelObject & model_object,
        MapObject & map_object,
        const PotentialAnalysisRequest & request);
    void SavePreparedModel(ModelObject & model_object, std::string_view saved_key_tag);

    void StudyAtomLocalFittingViaAlphaR(
        ModelObject & model_object,
        const std::vector<AtomObject *> & atom_list,
        const std::vector<double> & alpha_list,
        const std::filesystem::path & training_report_dir
    );
    void StudyAtomGroupFittingViaAlphaG(
        ModelObject & model_object,
        const std::vector<std::vector<AtomObject *>> & atom_list_set,
        const std::vector<double> & alpha_list,
        const std::filesystem::path & training_report_dir
    );

};

} // namespace rhbm_gem
