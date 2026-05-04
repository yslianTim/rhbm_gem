#pragma once

#include <Eigen/Dense>

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class RHBMTestCommand : public CommandBase<RHBMTestRequest>
{
public:
    RHBMTestCommand();
    ~RHBMTestCommand() = default;

private:
    void NormalizeAndValidateRequest() override;
    void ValidatePreparedRequest() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
