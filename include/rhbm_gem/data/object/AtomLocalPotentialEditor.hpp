#pragma once

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem {

class LocalPotentialEntry;
class ModelAnalysisEditor;

class AtomLocalPotentialEditor
{
    LocalPotentialEntry * m_entry;

public:
    void SetSamplingEntries(LocalPotentialSampleList value);
    void SetGaussianResult(LocalGaussianResult value);
    void SetAlphaR(double value);

private:
    explicit AtomLocalPotentialEditor(LocalPotentialEntry & entry);
    friend class ModelAnalysisEditor;

};

} // namespace rhbm_gem
