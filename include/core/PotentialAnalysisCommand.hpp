#pragma once

#include <iostream>
#include <string>
#include "CommandBase.hpp"

class PotentialAnalysisCommand : public CommandBase
{
    int m_thread_size, m_sampling_size;
    double m_sampling_range_min, m_sampling_range_max;
    double m_fit_range_min, m_fit_range_max;
    double m_alpha_r, m_alpha_g;
    std::string m_database_path, m_model_file_path, m_map_file_path;
    std::string m_pick_chain_id, m_pick_indicator;
    std::string m_pick_residue, m_pick_element, m_pick_remoteness, m_pick_branch;
    std::string m_veto_chain_id, m_veto_indicator;
    std::string m_veto_residue, m_veto_element, m_veto_remoteness, m_veto_branch;


public:
    PotentialAnalysisCommand(void);
    ~PotentialAnalysisCommand() = default;
    void Execute(void) override;

    void SetThreadSize(int value) { m_thread_size = value; }
    void SetSamplingSize(int value) { m_sampling_size = value; }
    void SetSamplingRangeMinimum(double value) { m_sampling_range_min = value; }
    void SetSamplingRangeMaximum(double value) { m_sampling_range_max = value; }
    void SetFitRangeMinimum(double value) { m_fit_range_min = value; }
    void SetFitRangeMaximum(double value) { m_fit_range_max = value; }
    void SetAlphaR(double value) { m_alpha_r = value; }
    void SetAlphaG(double value) { m_alpha_g = value; }
    void SetDatabasePath(const std::string & path) { m_database_path = path; }
    void SetModelFilePath(const std::string & path) { m_model_file_path = path; }
    void SetMapFilePath(const std::string & path) { m_map_file_path = path; }
    void SetPickChainID(const std::string & value) { m_pick_chain_id = value; }
    void SetPickIndicator(const std::string & value) { m_pick_indicator = value; }
    void SetPickResidueType(const std::string & value) { m_pick_residue = value; }
    void SetPickElementType(const std::string & value) { m_pick_element = value; }
    void SetPickRemotenessType(const std::string & value) { m_pick_remoteness = value; }
    void SetPickBranchType(const std::string & value) { m_pick_branch = value; }
    void SetVetoChainID(const std::string & value) { m_veto_chain_id = value; }
    void SetVetoIndicator(const std::string & value) { m_veto_indicator = value; }
    void SetVetoResidueType(const std::string & value) { m_veto_residue = value; }
    void SetVetoElementType(const std::string & value) { m_veto_element = value; }
    void SetVetoRemotenessType(const std::string & value) { m_veto_remoteness = value; }
    void SetVetoBranchType(const std::string & value) { m_veto_branch = value; }
};