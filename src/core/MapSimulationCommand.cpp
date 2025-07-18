#include "MapSimulationCommand.hpp"
#include "DataObjectManager.hpp"
#include "MapSimulationVisitor.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

namespace {
CommandRegistrar<MapSimulationCommand> registrar_map_simulation{
    "map_simulation",
    "Run map simulation command"};
}

MapSimulationCommand::MapSimulationCommand(void)
{
    SetBlurringWidthList(m_options.blurring_width_list);
}

void MapSimulationCommand::RegisterCLIOptions(CLI::App * cmd)
{
    cmd->add_option("-a,--model", m_options.model_file_path,
        "Model file path")->required();
    cmd->add_option("-n,--name", m_options.map_file_name,
        "File name for output map files")->default_val(m_options.map_file_name);
    cmd->add_option("--potential-model", m_options.potential_model_choice,
        "Atomic potential model option")->default_val(m_options.potential_model_choice);
    cmd->add_option("--charge", m_options.partial_charge_choice,
        "Partial charge table option")->default_val(m_options.partial_charge_choice);
    cmd->add_option("-c,--cut-off", m_options.cutoff_distance,
        "Cutoff distance")->default_val(m_options.cutoff_distance);
    cmd->add_option("-g,--grid-spacing", m_options.grid_spacing,
        "Grid spacing")->default_val(m_options.grid_spacing);
    cmd->add_option("--blurring-width", m_options.blurring_width_list,
        "Blurring width (list) setting")->default_val(m_options.blurring_width_list);
    RegisterCommandOptions(cmd);
}

bool MapSimulationCommand::Execute(void)
{
    Logger::Log(LogLevel::Info, "MapSimulationCommand::Execute() called.");
    try
    {
        SetBlurringWidthList(m_options.blurring_width_list);
        Logger::Log(LogLevel::Info, "Total number of blurring width sets to be simulated: "
                    + std::to_string(m_blurring_width_list.size()));

        auto data_manager{ std::make_unique<DataObjectManager>() };
        data_manager->ProcessFile(m_options.model_file_path, "model");

        auto analyzer{ std::make_unique<MapSimulationVisitor>() };
        analyzer->SetModelObjectKeyTag("model");
        analyzer->SetFolderPath(m_options.folder_path);
        analyzer->SetMapFileName(m_options.map_file_name);
        analyzer->SetThreadSize(static_cast<unsigned int>(m_options.thread_size));
        analyzer->SetPotentialModelChoice(m_options.potential_model_choice);
        analyzer->SetPartialChargeChoice(m_options.partial_charge_choice);
        analyzer->SetCutoffDistance(m_options.cutoff_distance);
        analyzer->SetGridSpacing(m_options.grid_spacing);
        analyzer->SetBlurringWidthList(m_blurring_width_list);

        data_manager->Accept(analyzer.get());
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error, e.what());
        return false;
    }
    return true;
}

void MapSimulationCommand::SetBlurringWidthList(const std::string & value)
{
    m_blurring_width_list.clear();
    for (const auto & token : StringHelper::SplitStringLineFromDelimiter(value, ','))
    {
        m_blurring_width_list.emplace_back(std::stod(token));
    }
}
