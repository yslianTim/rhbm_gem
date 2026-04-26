#pragma once

#include <Eigen/Dense>

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class RHBMTestCommand : public CommandWithRequest<RHBMTestRequest>
{
public:
    RHBMTestCommand();
    ~RHBMTestCommand() = default;

private:
    void NormalizeRequest() override;
    void ValidateOptions() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
