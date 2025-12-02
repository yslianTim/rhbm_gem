#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

class ModelObject;
class MapObject;
class AtomObject;

class PotentialAnalysisCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
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
        std::string saved_key_tag{"model"};
    };

private:
    Options m_options;
    std::string m_model_key_tag, m_map_key_tag;
    std::shared_ptr<MapObject> m_map_object;
    std::shared_ptr<ModelObject> m_model_object;

public:
    PotentialAnalysisCommand(void);
    ~PotentialAnalysisCommand();
    bool Execute(void) override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    const CommandOptions & GetOptions(void) const override { return m_options; }
    CommandOptions & GetOptions(void) override { return m_options; }

    void SetTrainingAlphaFlag(bool value);
    void SetAsymmetryFlag(bool value);
    void SetSimulationFlag(bool value);
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

private:
    bool BuildDataObject(void);
    void UpdateModelObjectForSimulation(ModelObject * model_object);
    void RunMapObjectPreprocessing(void);
    void RunModelObjectPreprocessing(void);
    void RunAtomMapValueSampling(void);
    void RunBondMapValueSampling(void);
    void RunAtomGroupClassification(void);
    void RunBondGroupClassification(void);
    void RunAtomAlphaTraining(void);
    std::vector<double> TrainAlphaR(const AtomObject * atom, const size_t group_size, const std::vector<double> & alpha_list);
    std::vector<double> TrainAlphaG(const std::vector<AtomObject *> & atom_list, const size_t group_size, const std::vector<double> & alpha_list);
    double TrainUniversalAlphaR(const std::vector<AtomObject *> & atom_list, const size_t group_size, const std::vector<double> & alpha_list);
    double TrainUniversalAlphaG(const std::vector<AtomObject *> & atom_list, const size_t group_size, const std::vector<double> & alpha_list);
    void RunLocalAtomFitting(void);
    void RunLocalBondFitting(void);
    void RunAtomPotentialFitting(void);
    void RunBondPotentialFitting(void);
    void SaveDataObject(void);

};
