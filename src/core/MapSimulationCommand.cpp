#include "MapSimulationCommand.hpp"
#include "DataObjectManager.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "MapFileWriter.hpp"
#include "ScopeTimer.hpp"
#include "FilePathHelper.hpp"
#include "ElectricPotential.hpp"
#include "KDTreeAlgorithm.hpp"
#include "AminoAcidInfoHelper.hpp"
#include "ArrayStats.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

#include <algorithm>
#include <limits>

namespace {
CommandRegistrar<MapSimulationCommand> registrar_map_simulation{
    "map_simulation",
    "Run map simulation command"};
}

MapSimulationCommand::MapSimulationCommand(void) :
    CommandBase(), m_options{}, m_selected_atom_list{}, m_atom_charge_map{},
    m_model_object{ nullptr },
    m_atom_range_minimum{
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max() },
    m_atom_range_maximum{
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest() }
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::MapSimulationCommand() called");
}

MapSimulationCommand::~MapSimulationCommand()
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::~MapSimulationCommand() called");
}

void MapSimulationCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::RegisterCLIOptionsExtend() called");
    cmd->add_option_function<std::string>("-a,--model",
        [&](const std::string & value) { SetModelFilePath(value); },
        "Model file path")->required();
    cmd->add_option_function<std::string>("-n,--name",
        [&](const std::string & value) { SetMapFileName(value); },
        "File name for output map files")->default_val(m_options.map_file_name);
    std::map<std::string, PotentialModel> model_map
    {
        {"0", PotentialModel::SINGLE_GAUS},      {"single", PotentialModel::SINGLE_GAUS},
        {"1", PotentialModel::FIVE_GAUS_CHARGE}, {"five",   PotentialModel::FIVE_GAUS_CHARGE},
        {"2", PotentialModel::SINGLE_GAUS_USER}, {"user",   PotentialModel::SINGLE_GAUS_USER}
    };
    cmd->add_option_function<PotentialModel>("--potential-model",
        [&](PotentialModel value) { SetPotentialModelChoice(value); },
        "Atomic potential model option")
        ->default_val(PotentialModel::FIVE_GAUS_CHARGE)
        ->transform(CLI::CheckedTransformer(model_map, CLI::ignore_case));
    std::map<std::string, PartialCharge> charge_map
    {
        {"0", PartialCharge::NEUTRAL}, {"neutral", PartialCharge::NEUTRAL},
        {"1", PartialCharge::PARTIAL}, {"partial", PartialCharge::PARTIAL},
        {"2", PartialCharge::AMBER},   {"amber",   PartialCharge::AMBER}
    };
    cmd->add_option_function<PartialCharge>("--charge",
        [&](PartialCharge value) { SetPartialChargeChoice(value); },
        "Partial charge table option")
        ->default_val(PartialCharge::PARTIAL)
        ->transform(CLI::CheckedTransformer(charge_map, CLI::ignore_case));
    cmd->add_option_function<double>("-c,--cut-off",
        [&](double value) { SetCutoffDistance(value); },
        "Cutoff distance")->default_val(m_options.cutoff_distance);
    cmd->add_option_function<double>("-g,--grid-spacing",
        [&](double value) { SetGridSpacing(value); },
        "Grid spacing")->default_val(m_options.grid_spacing);
    cmd->add_option_function<std::string>("--blurring-width",
        [&](const std::string & value) { SetBlurringWidthList(value); },
        "Blurring width (list) setting")->default_val(m_options.blurring_width_list);
}

bool MapSimulationCommand::Execute(void)
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::Execute() called");
    if (BuildDataObject() == false) return false;
    CalculateAtomRange();
    RunMapSimulation();
    return true;
}

void MapSimulationCommand::SetPotentialModelChoice(PotentialModel value)
{
    m_options.potential_model_choice = value;
}

void MapSimulationCommand::SetPartialChargeChoice(PartialCharge value)
{
    m_options.partial_charge_choice = value;
}

void MapSimulationCommand::SetCutoffDistance(double value)
{
    m_options.cutoff_distance = value;
    if (m_options.cutoff_distance <= 0.0)
    {
        Logger::Log(LogLevel::Warning,
            "Cutoff distance must be positive, reset to default 5.0");
        m_options.cutoff_distance = 5.0;
    }
}

void MapSimulationCommand::SetModelFilePath(const std::filesystem::path & value)
{
    m_options.model_file_path = value;
    if (!FilePathHelper::EnsureFileExists(m_options.model_file_path, "Model file"))
    {
        Logger::Log(LogLevel::Error,
            "Model file does not exist: " + m_options.model_file_path.string());
        m_valiate_options = false;
    }
}

void MapSimulationCommand::SetMapFileName(const std::string & value)
{
    m_options.map_file_name = value;
}

void MapSimulationCommand::SetGridSpacing(double value)
{
    m_options.grid_spacing = value;
    if (m_options.grid_spacing <= 0.0)
    {
        Logger::Log(LogLevel::Warning,
            "Grid spacing must be positive, reset to default 0.5");
        m_options.grid_spacing = 0.5;
    }
}

void MapSimulationCommand::SetBlurringWidthList(const std::string & value)
{
    m_options.blurring_width_list = StringHelper::ParseListOption<double>(value, ',');
    for (auto width : m_options.blurring_width_list)
    {
        if (width <= 0.0)
        {
            Logger::Log(LogLevel::Warning,
                "Blurring width must be positive, erase current setting : "
                + std::to_string(width));
            m_options.blurring_width_list.erase(
                std::remove(
                    m_options.blurring_width_list.begin(),
                    m_options.blurring_width_list.end(), width),
                m_options.blurring_width_list.end()
            );
        }
    }
}

bool MapSimulationCommand::BuildDataObject(void)
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::BuildDataObject() called");
    ScopeTimer timer("MapSimulationCommand::BuildDataObject");
    try
    {
        auto data_manager{ GetDataManagerPtr() };
        data_manager->ProcessFile(m_options.model_file_path, "model");
        m_model_object = data_manager->GetTypedDataObject<ModelObject>("model");
        BuildAtomList(m_model_object.get());
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "MapSimulationCommand::BuildDataObject : " + std::string(e.what()));
        return false;
    }
    return true;
}

void MapSimulationCommand::RunMapSimulation(void)
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::RunMapSimulation() called");
    ScopeTimer timer("MapSimulationCommand::RunMapSimulation");
    
    Logger::Log(LogLevel::Info,
        "Total number of blurring width sets to be simulated: "
        + std::to_string(m_options.blurring_width_list.size()));
    
    auto map_object{ CreateMapObject() };
    for (auto & blurring_width : m_options.blurring_width_list)
    {
        auto map_key_tag{
            m_model_object->GetPdbID() + "_bw" +
            StringHelper::ToStringWithPrecision<double>(blurring_width, 2)
        };
        PopulateMapValueArray(map_object.get(), blurring_width);
        std::string file_name{ m_options.map_file_name + "_" + map_key_tag + ".map" };
        std::filesystem::path output_file_name{ m_options.folder_path / file_name };
        MapFileWriter writer{ output_file_name.string(), map_object.get() };
        writer.Write();
    }
}

void MapSimulationCommand::BuildAtomList(ModelObject * model_object)
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::BuildAtomList() called");
    m_selected_atom_list.clear();
    m_selected_atom_list.reserve(model_object->GetNumberOfAtom());
    m_atom_charge_map.clear();
    for (auto & atom : model_object->GetComponentsList())
    {
        if (atom->IsUnknownAtom() == true)
        {
            Logger::Log(LogLevel::Warning,
                "Unknown atom found in the model object: Serial-ID = "
                + std::to_string(atom->GetSerialID()) +
                ", this atom will be ignored in the simulation.");
            continue;
        }
        m_selected_atom_list.emplace_back(atom.get());
        m_atom_charge_map.emplace(atom->GetSerialID(), CalculateAtomCharge(atom.get()));
    }

    Logger::Log(LogLevel::Info,
        "Number of selected atoms to be simulated = "
        + std::to_string(m_selected_atom_list.size()) +" / "
        + std::to_string(model_object->GetNumberOfAtom()) + " atoms.");
}

double MapSimulationCommand::CalculateAtomCharge(AtomObject * atom) const
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::CalculateAtomCharge() called");
    switch (m_options.partial_charge_choice)
    {
        case PartialCharge::NEUTRAL:
            return 0.0;
        case PartialCharge::PARTIAL:
            return AminoAcidInfoHelper::GetPartialCharge(
                atom->GetResidue(),
                atom->GetElement(),
                atom->GetRemoteness(),
                atom->GetBranch(),
                atom->GetStructure());
        case PartialCharge::AMBER:
            return AminoAcidInfoHelper::GetPartialCharge(
                atom->GetResidue(),
                atom->GetElement(),
                atom->GetRemoteness(),
                atom->GetBranch(),
                atom->GetStructure(), true);
        default:
            Logger::Log(LogLevel::Error,
                "Invalid partial charge choice: "
                + std::to_string(static_cast<int>(m_options.partial_charge_choice)));
            break;
    }
    return 0.0;
}

void MapSimulationCommand::CalculateAtomRange(void)
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::CalculateAtomRange() called");
    if (m_selected_atom_list.empty())
    {
        Logger::Log(LogLevel::Warning, "No atoms selected. Atom range cannot be calculated.");
        return;
    }

    for (const auto & atom : m_selected_atom_list)
    {
        const auto & atom_position{ atom->GetPositionRef() };
        m_atom_range_minimum[0] = std::min(m_atom_range_minimum[0], atom_position[0]);
        m_atom_range_minimum[1] = std::min(m_atom_range_minimum[1], atom_position[1]);
        m_atom_range_minimum[2] = std::min(m_atom_range_minimum[2], atom_position[2]);
        m_atom_range_maximum[0] = std::max(m_atom_range_maximum[0], atom_position[0]);
        m_atom_range_maximum[1] = std::max(m_atom_range_maximum[1], atom_position[1]);
        m_atom_range_maximum[2] = std::max(m_atom_range_maximum[2], atom_position[2]);
    }

    m_atom_range_minimum[0] -= static_cast<float>(m_options.cutoff_distance);
    m_atom_range_minimum[1] -= static_cast<float>(m_options.cutoff_distance);
    m_atom_range_minimum[2] -= static_cast<float>(m_options.cutoff_distance);
    m_atom_range_maximum[0] += static_cast<float>(m_options.cutoff_distance);
    m_atom_range_maximum[1] += static_cast<float>(m_options.cutoff_distance);
    m_atom_range_maximum[2] += static_cast<float>(m_options.cutoff_distance);
}

std::unique_ptr<MapObject> MapSimulationCommand::CreateMapObject(void)
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::CreateMapObject() called");
    ScopeTimer timer("MapSimulationCommand::CreateMapObject");
    std::array<float, 3> grid_spacing{
        static_cast<float>(m_options.grid_spacing),
        static_cast<float>(m_options.grid_spacing),
        static_cast<float>(m_options.grid_spacing)
    };

    auto origin{ CalculateOrigin(grid_spacing) };
    auto grid_size{ CalculateGridSize(grid_spacing) };
    auto map_object{ std::make_unique<MapObject>(grid_size, grid_spacing, origin) };
    map_object->SetThreadSize(static_cast<int>(m_options.thread_size));

    auto voxel_size{ map_object->GetMapValueArraySize() };
    auto map_value_array{ std::make_unique<float[]>(voxel_size) };
    std::fill_n(map_value_array.get(), voxel_size, 0.0f);
    map_object->SetMapValueArray(std::move(map_value_array));
    map_object->BuildKDTreeRoot();

    return map_object;
}

void MapSimulationCommand::PopulateMapValueArray(MapObject * map_object, double blurring_width)
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::PopulateMapValueArray() called");
    ScopeTimer timer("MapSimulationCommand::PopulateMapValueArray");
    Logger::Log(LogLevel::Info,
        " /- Start map value array production with blurring width = "+
        StringHelper::ToStringWithPrecision<double>(blurring_width, 2)
    );

    auto electric_potential{ std::make_unique<ElectricPotential>() };
    electric_potential->SetBlurringWidth(blurring_width);
    electric_potential->SetModelChoice(static_cast<int>(m_options.potential_model_choice));

    auto voxel_size{ map_object->GetMapValueArraySize() };
    auto map_value_array{ std::make_unique<float[]>(voxel_size) };
    std::fill_n(map_value_array.get(), voxel_size, 0.0f);

    auto atom_size{ m_selected_atom_list.size() };
    auto kd_tree_root{ map_object->GetKDTreeRoot() };
    size_t atom_count{ 0 };
    std::vector<GridNode*> in_range_grid_node_list;

    Logger::Log(LogLevel::Info,
        " /- Total number of atoms to be processed: "+ std::to_string(atom_size) + " atoms."
    );
#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(m_options.thread_size) private(in_range_grid_node_list)
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        auto atom{ m_selected_atom_list[i] };
        auto charge{ m_atom_charge_map.at(atom->GetSerialID()) };
        auto element{ atom->GetElement() };
        auto atom_position{ atom->GetPosition() };
        MapObject query_map_object({1, 1, 1}, {1.0, 1.0, 1.0}, atom_position);
        GridNode query_node(0, &query_map_object);

        KDTreeAlgorithm<GridNode>::RangeSearch(
            kd_tree_root, &query_node, m_options.cutoff_distance, in_range_grid_node_list);
        
        for (auto & grid_node : in_range_grid_node_list)
        {
            auto distance{
                ArrayStats<float>::ComputeNorm(atom_position, grid_node->GetPosition())
            };
            map_value_array[grid_node->GetGridIndex()] += static_cast<float>(
                electric_potential->GetPotentialValue(element, distance, charge)
            );
        }

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, atom_size);
        }
    }

    map_object->SetMapValueArray(std::move(map_value_array));
    map_object->Update();
    map_object->Display();
}

std::array<int, 3> MapSimulationCommand::CalculateGridSize(
    const std::array<float, 3> & grid_spacing) const
{
    auto selected_atom_size{ m_selected_atom_list.size() };
    if (selected_atom_size == 0)
    {
        Logger::Log(LogLevel::Warning, "No atoms selected. Grid size is set to [1,1,1].");
        return { 1, 1, 1 };
    }

    return {
        static_cast<int>(
            std::ceil((m_atom_range_maximum[0] - m_atom_range_minimum[0]) / grid_spacing[0])),
        static_cast<int>(
            std::ceil((m_atom_range_maximum[1] - m_atom_range_minimum[1]) / grid_spacing[1])),
        static_cast<int>(
            std::ceil((m_atom_range_maximum[2] - m_atom_range_minimum[2]) / grid_spacing[2]))
    };
}

std::array<float, 3> MapSimulationCommand::CalculateOrigin(
    const std::array<float, 3> & grid_spacing) const
{
    Logger::Log(LogLevel::Debug, "MapSimulationCommand::CalculateOrigin() called");
    return {
        std::floor(m_atom_range_minimum[0] / grid_spacing[0]) * grid_spacing[0],
        std::floor(m_atom_range_minimum[1] / grid_spacing[1]) * grid_spacing[1],
        std::floor(m_atom_range_minimum[2] / grid_spacing[2]) * grid_spacing[2]
    };
}
