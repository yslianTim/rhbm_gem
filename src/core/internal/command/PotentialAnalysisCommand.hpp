#pragma once

#include <memory>
#include <string>
#include <filesystem>

#include "CommandBase.hpp"

namespace CLI
{
    class App;
}

namespace rhbm_gem {

class ModelObject;
class MapObject;
class AtomObject;
struct PotentialAnalysisRequest;
void BindPotentialAnalysisRequestOptions(CLI::App * command, PotentialAnalysisRequest & request);

struct PotentialAnalysisCommandOptions
{
    bool use_training_alpha{ false };
    bool is_asymmetry{ false };
    bool is_simulation{ false };
    int sampling_size{ 1500 };
    double sampling_range_min{ 0.0 };
    double sampling_range_max{ 1.5 };
    double sampling_height{ 0.1 };
    double fit_range_min{ 0.0 };
    double fit_range_max{ 1.0 };
    double alpha_r{ 0.1 };
    double alpha_g{ 0.2 };
    double resolution_simulation{ 0.0 };
    std::filesystem::path model_file_path;
    std::filesystem::path map_file_path;
    std::filesystem::path training_report_dir{};
    std::string saved_key_tag{"model"};
};

class PotentialAnalysisCommand : public CommandBase
{
private:
    PotentialAnalysisCommandOptions m_options{};
    std::string m_model_key_tag, m_map_key_tag;
    std::shared_ptr<MapObject> m_map_object;
    std::shared_ptr<ModelObject> m_model_object;

public:
    explicit PotentialAnalysisCommand(CommonOptionProfile profile);
    ~PotentialAnalysisCommand() override = default;
    void ApplyRequest(const PotentialAnalysisRequest & request);

private:
    void SetSimulatedMapResolution(double value);
    void SetFitRangeMinimum(double value);
    void SetFitRangeMaximum(double value);
    void SetAlphaR(double value);
    void SetAlphaG(double value);
    void SetModelFilePath(const std::filesystem::path & path);
    void SetMapFilePath(const std::filesystem::path & path);
    void SetSavedKeyTag(const std::string & tag);
    void SetSamplingSize(int value);
    void SetSamplingRangeMinimum(double value);
    void SetSamplingRangeMaximum(double value);
    void SetSamplingHeight(double value);
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
