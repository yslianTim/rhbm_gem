#include "MapSimulationCommand.hpp"
#include "DataObjectManager.hpp"
#include "MapSimulationVisitor.hpp"
#include "Logger.hpp"

#include <sstream>

MapSimulationCommand::MapSimulationCommand(void) :
    m_thread_size{ 1 }, m_potential_model_choice{ 1 },
    m_partial_charge_choice{ 1 }, m_cutoff_distance{ 5.0 },
    m_grid_spacing{ 0.5 }, m_model_file_path{ "" },
    m_folder_path{ "" }, m_map_file_name{ "sim_map" },
    m_blurring_width_list{ 0.50 }
{

}

MapSimulationCommand::MapSimulationCommand(
    const Options & options, const GlobalOptions & globals) :
    m_thread_size{ globals.thread_size },
    m_potential_model_choice{ options.potential_model_choice },
    m_partial_charge_choice{ options.partial_charge_choice },
    m_cutoff_distance{ options.cutoff_distance },
    m_grid_spacing{ options.grid_spacing },
    m_model_file_path{ options.model_file_path },
    m_folder_path{ globals.folder_path },
    m_map_file_name{ options.map_file_name }
{
    SetBlurringWidthList(options.blurring_width_list);
}

void MapSimulationCommand::Execute(void)
{
    Logger::Log(LogLevel::Info, "MapSimulationCommand::Execute() called.");
    Logger::Log(LogLevel::Info, "Total number of blurring width sets to be simulated: "
                + std::to_string(m_blurring_width_list.size()));

    auto data_manager{ std::make_unique<DataObjectManager>() };
    data_manager->ProcessFile(m_model_file_path, "model");

    auto analyzer{ std::make_unique<MapSimulationVisitor>() };
    analyzer->SetModelObjectKeyTag("model");
    analyzer->SetFolderPath(m_folder_path);
    analyzer->SetMapFileName(m_map_file_name);
    analyzer->SetThreadSize(static_cast<unsigned int>(m_thread_size));
    analyzer->SetPotentialModelChoice(m_potential_model_choice);
    analyzer->SetPartialChargeChoice(m_partial_charge_choice);
    analyzer->SetCutoffDistance(m_cutoff_distance);
    analyzer->SetGridSpacing(m_grid_spacing);
    analyzer->SetBlurringWidthList(m_blurring_width_list);

    data_manager->Accept(analyzer.get());
}

void MapSimulationCommand::SetBlurringWidthList(const std::string & value)
{
    m_blurring_width_list.clear();
    std::stringstream ss(value);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        if (segment == "") continue;
        m_blurring_width_list.emplace_back(std::stod(segment));
    }
}
