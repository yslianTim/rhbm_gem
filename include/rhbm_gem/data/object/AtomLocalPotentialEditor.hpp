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
    void SetRawSamplingEntries(LocalPotentialSampleList value);
    void SetPeelingSamplingEntries(LocalPotentialSampleList value);
    void SetNeighborCountForPeeling(int value);
    void SetGaussianResult(LocalFittingStage stage, LocalGaussianResult value);
    void SetAlphaR(LocalFittingStage stage, double value);

private:
    explicit AtomLocalPotentialEditor(LocalPotentialEntry & entry);
    friend class ModelAnalysisEditor;

};

} // namespace rhbm_gem
