#pragma once

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class ResultDumpCommand : public CommandBase<ResultDumpRequest>
{
public:
    ResultDumpCommand();
    ~ResultDumpCommand() override = default;

private:
    void NormalizeAndValidateRequest() override;
    void ValidatePreparedRequest() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
