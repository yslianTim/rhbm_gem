#pragma once

#include <memory>

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class ModelObject;
class MapObject;
class AtomObject;

class PotentialAnalysisCommand : public CommandWithRequest<PotentialAnalysisRequest>
{
private:
    struct AtomSamplingOptions;

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
    void UpdateModelObjectForSimulation(
        ModelObject * model_object,
        double simulated_map_resolution);
    void RunMapObjectPreprocessing();
    void RunModelObjectPreprocessing(bool asymmetry_flag);
    void RunAtomSamplingWorkflow(const AtomSamplingOptions & options);
    void RunAtomGroupingWorkflow();
    void RunAtomAlphaTraining(const std::filesystem::path & training_report_dir);
    double TrainUniversalAlphaR(
        const std::vector<AtomObject *> & atom_list,
        const size_t subset_size,
        const std::vector<double> & alpha_list
    );
    double TrainUniversalAlphaG(
        const std::vector<std::vector<AtomObject *>> & atom_list_set,
        const size_t subset_size,
        const std::vector<double> & alpha_list
    );
    void RunLocalAtomFittingWorkflow(double universal_alpha_r);
    void RunAtomPotentialFitting(bool training_alpha_flag, double fallback_alpha_g);
    void RunExperimentalBondWorkflowIfEnabled(const PotentialAnalysisRequest & request);
    void SavePreparedModel(std::string_view saved_key_tag);

    void StudyAtomLocalFittingViaAlphaR(
        const std::vector<AtomObject *> & atom_list,
        const std::vector<double> & alpha_list,
        const std::filesystem::path & training_report_dir
    );
    void StudyAtomGroupFittingViaAlphaG(
        const std::vector<std::vector<AtomObject *>> & atom_list_set,
        const std::vector<double> & alpha_list,
        const std::filesystem::path & training_report_dir
    );

};

} // namespace rhbm_gem
