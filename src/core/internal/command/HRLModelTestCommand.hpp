#pragma once

#include <Eigen/Dense>

#include "CommandBase.hpp"

namespace rhbm_gem {

class HRLModelTestCommand : public CommandWithRequest<HRLModelTestRequest>
{
public:
    explicit HRLModelTestCommand(CommonOptionProfile profile);
    ~HRLModelTestCommand() = default;

private:
    void NormalizeRequest() override;
    void ValidateOptions() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
