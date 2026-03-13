#pragma once

#include "command/HRLModelTestCommand.hpp"

namespace rhbm_gem::detail {

struct HRLModelTestWorkflowContext
{
    const HRLModelTestCommandOptions & options;
};

bool RunHRLModelTestWorkflow(const HRLModelTestWorkflowContext & context);

} // namespace rhbm_gem::detail
