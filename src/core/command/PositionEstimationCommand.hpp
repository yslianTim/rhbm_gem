#pragma once

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class PositionEstimationCommand : public CommandBase<PositionEstimationRequest>
{
public:
    PositionEstimationCommand();
    ~PositionEstimationCommand() override = default;

private:
    void NormalizeAndValidateRequest() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
