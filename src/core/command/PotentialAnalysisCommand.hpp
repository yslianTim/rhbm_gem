#pragma once

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class PotentialAnalysisCommand : public CommandBase<PotentialAnalysisRequest>
{
public:
    PotentialAnalysisCommand();
    ~PotentialAnalysisCommand() override = default;

private:
    void NormalizeAndValidateRequest() override;
    void ValidatePreparedRequest() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
