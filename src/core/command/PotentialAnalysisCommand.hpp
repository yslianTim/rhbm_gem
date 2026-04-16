#pragma once

#include <memory>

#include <Eigen/Dense>

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class ModelObject;
class MapObject;
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
        const PotentialAnalysisRequest & request);
    void RunLocalPotentialFitting(ModelObject & model_object, double alpha_r);
    void RunAtomPotentialFittingWorkflow(
        ModelObject & model_object,
        const PotentialAnalysisRequest & request);
    void SavePreparedModel(ModelObject & model_object, std::string_view saved_key_tag);

};

} // namespace rhbm_gem
