#pragma once

#include "detail/CommandBase.hpp"

namespace rhbm_gem {

class MapVisualizationCommand : public CommandBase<MapVisualizationRequest>
{
public:
    MapVisualizationCommand();
    ~MapVisualizationCommand() override = default;

private:
    void NormalizeAndValidateRequest() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
