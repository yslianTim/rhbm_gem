#pragma once

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class MapSimulationCommand : public CommandBase<MapSimulationRequest>
{
public:
    MapSimulationCommand();
    ~MapSimulationCommand() override = default;

private:
    void NormalizeAndValidateRequest() override;
    void ValidatePreparedRequest() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
