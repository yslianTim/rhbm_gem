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
    bool BuildDataObject();
    void UpdateModelObjectForSimulation(ModelObject * model_object);
    void RunMapObjectPreprocessing();
    void RunModelObjectPreprocessing();
    void RunAtomMapValueSampling();
    void RunAtomGroupClassification();
    void RunAtomAlphaTraining();
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
    void RunLocalAtomFitting(double universal_alpha_r);
    void RunAtomPotentialFitting();
    void RunExperimentalBondWorkflowIfEnabled();
    void SaveDataObject();

    void StudyAtomLocalFittingViaAlphaR(
        const std::vector<AtomObject *> & atom_list,
        const std::vector<double> & alpha_list
    );
    void StudyAtomGroupFittingViaAlphaG(
        const std::vector<std::vector<AtomObject *>> & atom_list_set,
        const std::vector<double> & alpha_list
    );

};

} // namespace rhbm_gem
