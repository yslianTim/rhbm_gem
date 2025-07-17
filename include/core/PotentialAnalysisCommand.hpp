#pragma once

#include <memory>
#include <string>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "SphereSampler.hpp"

class ModelObject;
class PotentialAnalysisCommand : public CommandBase
{
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

private:
    Options m_options{};
    std::unique_ptr<SphereSampler> m_sphere_sampler;

public:
    PotentialAnalysisCommand(void);
    ~PotentialAnalysisCommand() = default;
    void Execute(void) override;
    void RegisterCLIOptions(CLI::App * cmd) override;
    void SetGlobalOptions(const GlobalOptions & globals) override;

    void SetAsymmetryFlag(bool value) { m_options.is_asymmetry = value; }
    void SetSimulationFlag(bool value) { m_options.is_simulation = value; }
    void SetSimulatedMapResolution(double value) { m_options.resolution_simulation = value; }
    void SetFitRangeMinimum(double value) { m_options.fit_range_min = value; }
    void SetFitRangeMaximum(double value) { m_options.fit_range_max = value; }
    void SetAlphaR(double value) { m_options.alpha_r = value; }
    void SetAlphaG(double value) { m_options.alpha_g = value; }
    void SetDatabasePath(const std::string & path) { m_globals.database_path = path; }
    void SetModelFilePath(const std::string & path) { m_options.model_file_path = path; }
    void SetMapFilePath(const std::string & path) { m_options.map_file_path = path; }
    void SetSavedKeyTag(const std::string & tag) { m_options.saved_key_tag = tag; }
    void SetThreadSize(int value);
    void SetSamplingSize(int value);
    void SetSamplingRangeMinimum(double value);
    void SetSamplingRangeMaximum(double value);

private:
    void UpdateModelObjectForSimulation(ModelObject * model_object);

};
