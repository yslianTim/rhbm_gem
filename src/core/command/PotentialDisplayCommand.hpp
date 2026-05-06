#pragma once

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class PotentialDisplayCommand : public CommandBase<PotentialDisplayRequest>
{
public:
    PotentialDisplayCommand();
    ~PotentialDisplayCommand() override = default;

private:
    void NormalizeAndValidateRequest() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
