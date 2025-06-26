#include "MapSimulationCommand.hpp"
#include "DataObjectManager.hpp"
#include "MapSimulationVisitor.hpp"

#include <iostream>
#include <sstream>

MapSimulationCommand::MapSimulationCommand(void)
{

}

MapSimulationCommand::~MapSimulationCommand()
{

}

void MapSimulationCommand::Execute(void)
{
    std::cout << "MapSimulationCommand::Execute() called." << std::endl;
    std::cout << "Total number of blurring width sets to be simulated: "
              << m_blurring_width_list.size() << std::endl;

    auto data_manager{ std::make_unique<DataObjectManager>() };
    data_manager->ProcessFile(m_model_file_path, "model");

    auto analyzer{ std::make_unique<MapSimulationVisitor>() };
    analyzer->SetModelObjectKeyTag("model");
    analyzer->SetFolderPath(m_folder_path);
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
    std::stringstream ss(value);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        if (segment == "") continue;
        m_blurring_width_list.emplace_back(std::stod(segment));
    }
}
