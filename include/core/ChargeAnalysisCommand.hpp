#pragma once

#include <memory>
#include <string>
#include "CommandBase.hpp"

class SphereSampler;
class ChargeAnalysisCommand : public CommandBase
{
    int m_thread_size;
    bool m_is_asymmetry;
    double m_fit_range_min, m_fit_range_max;
    double m_alpha_r, m_alpha_g;
    std::string m_database_path, m_model_file_path, m_map_file_path;
    std::string m_sim_neutral_map_file_path, m_sim_positive_map_file_path, m_sim_negative_map_file_path;
    std::string m_saved_key_tag;
    std::unique_ptr<SphereSampler> m_sphere_sampler;

public:
    ChargeAnalysisCommand(void);
    ~ChargeAnalysisCommand();
    void Execute(void) override;

    void SetAsymmetryFlag(bool value) { m_is_asymmetry = value; }
    void SetFitRangeMinimum(double value) { m_fit_range_min = value; }
    void SetFitRangeMaximum(double value) { m_fit_range_max = value; }
    void SetAlphaR(double value) { m_alpha_r = value; }
    void SetAlphaG(double value) { m_alpha_g = value; }
    void SetDatabasePath(const std::string & path) { m_database_path = path; }
    void SetModelFilePath(const std::string & path) { m_model_file_path = path; }
    void SetMapFilePath(const std::string & path) { m_map_file_path = path; }
    void SetSimulatedNeutralMapFilePath(const std::string & path) { m_sim_neutral_map_file_path = path; }
    void SetSimulatedPositiveMapFilePath(const std::string & path) { m_sim_positive_map_file_path = path; }
    void SetSimulatedNegativeMapFilePath(const std::string & path) { m_sim_negative_map_file_path = path; }
    void SetSavedKeyTag(const std::string & tag) { m_saved_key_tag = tag; }
    void SetThreadSize(int value);
    void SetSamplingSize(int value);
    void SetSamplingRangeMinimum(double value);
    void SetSamplingRangeMaximum(double value);

};