#pragma once

#include "command/PotentialAnalysisCommand.hpp"
#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>

namespace rhbm_gem::detail {

inline HRLExecutionOptions MakePotentialAnalysisExecutionOptions(
    const PotentialAnalysisCommandOptions & options,
    bool quiet_mode)
{
    HRLExecutionOptions execution_options;
    execution_options.quiet_mode = quiet_mode;
    execution_options.thread_size = options.thread_size;
    return execution_options;
}

} // namespace rhbm_gem::detail
