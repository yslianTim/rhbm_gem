#pragma once

namespace rhbm_gem {

enum class CommandId
{
    PotentialAnalysis = 0,
    PotentialDisplay,
    ResultDump,
    MapSimulation,
    MapVisualization,
    PositionEstimation,
    ModelTest
};

enum class DatabaseUsage
{
    Required = 0,
    NotUsed
};

enum class BindingExposure
{
    PythonPublic = 0,
    CliOnly
};

} // namespace rhbm_gem
