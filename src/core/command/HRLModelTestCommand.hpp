#pragma once

#include <Eigen/Dense>

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class HRLModelTestCommand : public CommandWithRequest<HRLModelTestRequest>
{
public:
    HRLModelTestCommand();
    ~HRLModelTestCommand() = default;

private:
    void NormalizeRequest() override;
    void ValidateOptions() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
