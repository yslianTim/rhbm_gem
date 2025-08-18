#include "PositionEstimationCommand.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "KDTreeAlgorithm.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

#include <random>
#include <memory>
#include <vector>

namespace {
CommandRegistrar<PositionEstimationCommand> registrar_model_test{
    "position_estimation",
    "Run atom position estimation"};
}

PositionEstimationCommand::PositionEstimationCommand(void) :
    CommandBase(), m_options{}
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::PositionEstimationCommand() called.");
}

void PositionEstimationCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::RegisterCLIOptionsExtend() called.");
    cmd->add_option("-m,--map", m_options.map_file_path,
        "Map file path")->required();
    cmd->add_option("--alpha", m_options.alpha,
        "Alpha value for robust regression")->default_val(m_options.alpha);
    cmd->add_option("--iter", m_options.iteration_count,
        "Iteration count for estimation")->default_val(m_options.iteration_count);
    cmd->add_option("-k,--knn", m_options.knn_size,
        "KNN size for estimation")->default_val(m_options.knn_size);
}

bool PositionEstimationCommand::Execute(void)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::Execute() called.");
    auto data_manager{ GetDataManagerPtr() };
    data_manager->SetDatabaseManager(m_options.database_path);
    try
    {
        data_manager->ProcessFile(m_options.map_file_path, "map");
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PositionEstimationCommand::Execute() : " + std::string(e.what()));
        return false;
    }

    auto map_object{ data_manager->GetTypedDataObject<MapObject>("map") };
    map_object->MapValueArrayNormalization();

    RunMapValueConvergence(map_object.get());

    //data_manager->SaveDataObject("map", m_options.saved_key_tag);
    return true;
}

bool PositionEstimationCommand::ValidateOptions(void) const
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::ValidateOptions() called");

    return true;
}

void PositionEstimationCommand::RunMapValueConvergence(MapObject * map_object)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::RunMapValueConvergence() called");
    map_object->Display();

    for (int i = 0; i < m_options.iteration_count; i++)
    {
        Logger::Log(LogLevel::Debug, "Iteration: " + std::to_string(i + 1));
    }
}
