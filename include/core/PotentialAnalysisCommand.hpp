#pragma once

#include <memory>
#include <string>
#include "CommandBase.hpp"
#include "GlobalOptions.hpp"

class SphereSampler;
class ModelObject;
class PotentialAnalysisCommand : public CommandBase
{
    int m_thread_size;
    bool m_is_asymmetry, m_is_simulation;
    double m_fit_range_min, m_fit_range_max;
    double m_alpha_r, m_alpha_g;
    double m_simulated_map_resolution;
    std::string m_database_path, m_model_file_path, m_map_file_path;
    std::string m_saved_key_tag;
    std::unique_ptr<SphereSampler> m_sphere_sampler;

public:
    struct Options
    {
        bool is_asymmetry{ false };
        bool is_simulation{ false };
        int sampling_size{ 1500 };
        double sampling_range_min{ 0.0 };
        double sampling_range_max{ 1.5 };
        double fit_range_min{ 0.0 };
        double fit_range_max{ 1.0 };
        double alpha_r{ 0.1 };
        double alpha_g{ 0.2 };
        double resolution_simulation{ 0.0 };
        std::string model_file_path;
        std::string map_file_path;
        std::string saved_key_tag{"model"};
    };
    PotentialAnalysisCommand(void);
    PotentialAnalysisCommand(const Options & options, const GlobalOptions & globals);
    ~PotentialAnalysisCommand();
    void Execute(void) override;

    void SetAsymmetryFlag(bool value) { m_is_asymmetry = value; }
    void SetSimulationFlag(bool value) { m_is_simulation = value; }
    void SetSimulatedMapResolution(double value) { m_simulated_map_resolution = value; }
    void SetFitRangeMinimum(double value) { m_fit_range_min = value; }
    void SetFitRangeMaximum(double value) { m_fit_range_max = value; }
    void SetAlphaR(double value) { m_alpha_r = value; }
    void SetAlphaG(double value) { m_alpha_g = value; }
    void SetDatabasePath(const std::string & path) { m_database_path = path; }
    void SetModelFilePath(const std::string & path) { m_model_file_path = path; }
    void SetMapFilePath(const std::string & path) { m_map_file_path = path; }
    void SetSavedKeyTag(const std::string & tag) { m_saved_key_tag = tag; }
    void SetThreadSize(int value);
    void SetSamplingSize(int value);
    void SetSamplingRangeMinimum(double value);
    void SetSamplingRangeMaximum(double value);

private:
    void UpdateModelObjectForSimulation(ModelObject * model_object);

};
