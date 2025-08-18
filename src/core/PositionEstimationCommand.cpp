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
    CommandBase(), m_options{}, m_kd_tree_root{ nullptr }
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

    BuildVoxelList(map_object.get());
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

void PositionEstimationCommand::BuildVoxelList(MapObject * map_object)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::BuildVoxelList() called");
    ScopeTimer timer("PositionEstimationCommand::BuildVoxelList");
    m_selected_voxel_list.clear();
    size_t positive_count{
        static_cast<size_t>(std::count_if(
            map_object->GetMapValueArray(),
            map_object->GetMapValueArray() + map_object->GetMapValueArraySize(),
            [](float v) { return v > 0.0f; })
        )
    };
    m_selected_voxel_list.reserve(positive_count);

    for (size_t i = 0; i < map_object->GetMapValueArraySize(); i++)
    {
        auto position{ map_object->GetGridPosition(i) };
        float value{ map_object->GetMapValue(i) };
        if (value <= 0.0f) continue; // Skip negative values
        m_selected_voxel_list.emplace_back(position, value);
    }
    m_selected_voxel_list.shrink_to_fit();

    Logger::Log(LogLevel::Info,
        "Number of selected voxels to be estimated = "
        + std::to_string(m_selected_voxel_list.size()) + " / "
        + std::to_string(map_object->GetMapValueArraySize()) + " voxels."
    );

    std::vector<VoxelNode *> voxel_ptrs;
    voxel_ptrs.reserve(m_selected_voxel_list.size());
    for (auto & voxel : m_selected_voxel_list)
    {
        voxel_ptrs.emplace_back(&voxel);
    }

    m_kd_tree_root = KDTreeAlgorithm<VoxelNode>::BuildKDTree(voxel_ptrs, 0);
}
