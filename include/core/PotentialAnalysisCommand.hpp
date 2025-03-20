#pragma once

#include <iostream>
#include <memory>
#include <string>
#include "CommandBase.hpp"

class AtomSelector;
class SphereSampler;
class PotentialAnalysisCommand : public CommandBase
{
    int m_thread_size;
    double m_fit_range_min, m_fit_range_max;
    double m_alpha_r, m_alpha_g;
    std::string m_database_path, m_model_file_path, m_map_file_path;
    std::string m_saved_key_tag;
    std::shared_ptr<AtomSelector> m_atom_selector;
    std::shared_ptr<SphereSampler> m_sphere_sampler;

public:
    PotentialAnalysisCommand(void);
    ~PotentialAnalysisCommand();
    void Execute(void) override;

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
    void SetPickChainID(const std::string & value);
    void SetPickIndicator(const std::string & value);
    void SetPickResidueType(const std::string & value);
    void SetPickElementType(const std::string & value);
    void SetPickRemotenessType(const std::string & value);
    void SetPickBranchType(const std::string & value);
    void SetVetoChainID(const std::string & value);
    void SetVetoIndicator(const std::string & value);
    void SetVetoResidueType(const std::string & value);
    void SetVetoElementType(const std::string & value);
    void SetVetoRemotenessType(const std::string & value);
    void SetVetoBranchType(const std::string & value);
};