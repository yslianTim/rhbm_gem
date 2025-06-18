#pragma once

#include <memory>
#include <string>
#include "CommandBase.hpp"

class ChargeAnalysisCommand : public CommandBase
{
    int m_thread_size;
    double m_fit_range_min, m_fit_range_max;
    double m_alpha_r, m_alpha_g;
    std::string m_database_path;
    std::string m_sim_neutral_model_key_tag, m_sim_positive_model_key_tag, m_sim_negative_model_key_tag;
    std::string m_model_key_tag;

public:
    ChargeAnalysisCommand(void);
    ~ChargeAnalysisCommand();
    void Execute(void) override;

    void SetFitRangeMinimum(double value) { m_fit_range_min = value; }
    void SetFitRangeMaximum(double value) { m_fit_range_max = value; }
    void SetAlphaR(double value) { m_alpha_r = value; }
    void SetAlphaG(double value) { m_alpha_g = value; }
    void SetDatabasePath(const std::string & path) { m_database_path = path; }
    void SetSimulatedNeutralModelKeyTag(const std::string & path) { m_sim_neutral_model_key_tag = path; }
    void SetSimulatedPositiveModelKeyTag(const std::string & path) { m_sim_positive_model_key_tag = path; }
    void SetSimulatedNegativeModelKeyTag(const std::string & path) { m_sim_negative_model_key_tag = path; }
    void SetModelKeyTag(const std::string & value) { m_model_key_tag = value; }
    void SetThreadSize(int value) { m_thread_size = value; }

};
