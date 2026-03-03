#include "MapSimulationCommand.hpp"
#include "CommandOptionBinding.hpp"
#include "DataObjectManager.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "MapFileWriter.hpp"
#include "ScopeTimer.hpp"
#include "FilePathHelper.hpp"
#include "ElectricPotential.hpp"
#include "KDTreeAlgorithm.hpp"
#include "ComponentHelper.hpp"
#include "ArrayStats.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"
#include "OptionEnumTraits.hpp"

#include <algorithm>
#include <limits>

namespace {
constexpr std::string_view kModelKey{ "model" };

std::string SerializeBlurringWidths(const std::vector<double> & widths)
{
    std::string output;
    for (size_t index = 0; index < widths.size(); ++index)
    {
        if (index != 0) output += ",";
        output += StringHelper::ToStringWithPrecision<double>(widths[index], 2);
    }
    return output;
}
}

namespace rhbm_gem {

MapSimulationCommand::MapSimulationCommand() :
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
}

MapSimulationCommand::~MapSimulationCommand()
{
}

void MapSimulationCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    command_cli::AddPathOption(
        cmd, "-a,--model",
        [&](const std::filesystem::path & value) { SetModelFilePath(value); },
        "Model file path",
        std::nullopt,
        true);
    command_cli::AddStringOption(
        cmd, "-n,--name",
        [&](const std::string & value) { SetMapFileName(value); },
        "File name for output map files",
        m_options.map_file_name);
    command_cli::AddEnumOption<PotentialModel>(
        cmd, "--potential-model",
        [&](PotentialModel value) { SetPotentialModelChoice(value); },
        "Atomic potential model option",
        PotentialModel::FIVE_GAUS_CHARGE);
    command_cli::AddEnumOption<PartialCharge>(
        cmd, "--charge",
        [&](PartialCharge value) { SetPartialChargeChoice(value); },
        "Partial charge table option",
        PartialCharge::PARTIAL);
    command_cli::AddScalarOption<double>(
        cmd, "-c,--cut-off",
        [&](double value) { SetCutoffDistance(value); },
        "Cutoff distance",
        m_options.cutoff_distance);
    command_cli::AddScalarOption<double>(
        cmd, "-g,--grid-spacing",
        [&](double value) { SetGridSpacing(value); },
        "Grid spacing",
        m_options.grid_spacing);
    command_cli::AddStringOption(
        cmd, "--blurring-width",
        [&](const std::string & value) { SetBlurringWidthList(value); },
        "Blurring width (list) setting",
        SerializeBlurringWidths(m_options.blurring_width_list));
}

bool MapSimulationCommand::ExecuteImpl()
{
    if (BuildDataObject() == false) return false;
    CalculateAtomRange();
    RunMapSimulation();
    return true;
}

void MapSimulationCommand::ValidateOptions()
{
    ResetPrepareIssues("--blurring-width");
    if (m_options.blurring_width_list.empty())
    {
        AddValidationError("--blurring-width",
            "At least one positive blurring width is required.");
    }
}

void MapSimulationCommand::ResetRuntimeState()
{
    m_selected_atom_list.clear();
    m_atom_charge_map.clear();
    m_model_object.reset();
    m_atom_range_minimum = {
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max()
    };
    m_atom_range_maximum = {
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest()
    };
}

void MapSimulationCommand::SetPotentialModelChoice(PotentialModel value)
{
    MutateOptions([&]() { m_options.potential_model_choice = value; });
}

void MapSimulationCommand::SetPartialChargeChoice(PartialCharge value)
{
    MutateOptions([&]() { m_options.partial_charge_choice = value; });
}

void MapSimulationCommand::SetCutoffDistance(double value)
{
    SetNormalizedScalarOption(
        m_options.cutoff_distance,
        value,
        "--cut-off",
        [](double candidate) { return candidate > 0.0; },
        5.0,
        "Cutoff distance must be positive, reset to default 5.0");
}

void MapSimulationCommand::SetModelFilePath(const std::filesystem::path & value)
{
    SetRequiredExistingPathOption(m_options.model_file_path, value, "--model", "Model file");
}

void MapSimulationCommand::SetMapFileName(const std::string & value)
{
    MutateOptions([&]() { m_options.map_file_name = value; });
}

void MapSimulationCommand::SetGridSpacing(double value)
{
    SetNormalizedScalarOption(
        m_options.grid_spacing,
        value,
        "--grid-spacing",
        [](double candidate) { return candidate > 0.0; },
        0.5,
        "Grid spacing must be positive, reset to default 0.5");
}

void MapSimulationCommand::SetBlurringWidthList(const std::string & value)
{
    MutateOptions([&]()
    {
        ResetParseIssues("--blurring-width");
        const auto parsed_list{ StringHelper::ParseListOption<double>(value, ',') };
        m_options.blurring_width_list.clear();
        m_options.blurring_width_list.reserve(parsed_list.size());
        for (const auto width : parsed_list)
        {
            if (width <= 0.0)
            {
                AddNormalizationWarning("--blurring-width",
                    "Blurring width must be positive, dropping current setting: "
                        + std::to_string(width));
                continue;
            }
            m_options.blurring_width_list.push_back(width);
        }
    });
}

bool MapSimulationCommand::BuildDataObject()
{
    ScopeTimer timer("MapSimulationCommand::BuildDataObject");
    try
    {
        m_model_object = ProcessTypedFile<ModelObject>(
            m_options.model_file_path, kModelKey, "model file");
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

void MapSimulationCommand::RunMapSimulation()
{
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
        const auto output_file_name{
            BuildOutputPath(m_options.map_file_name + "_" + map_key_tag, ".map")
        };
        MapFileWriter writer{ output_file_name.string(), map_object.get() };
        writer.Write();
    }
}

void MapSimulationCommand::BuildAtomList(ModelObject * model_object)
{
    m_selected_atom_list.clear();
    m_selected_atom_list.reserve(model_object->GetNumberOfAtom());
    m_atom_charge_map.clear();
    for (auto & atom : model_object->GetAtomList())
    {
        if (atom->IsUnknownAtom() == true)
        {
            //Logger::Log(LogLevel::Warning,
            //    "Unknown atom found in the model object: Serial-ID = "
            //    + std::to_string(atom->GetSerialID()) +
            //    ", this atom will be ignored in the simulation.");
            //continue;
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
    switch (m_options.partial_charge_choice)
    {
        case PartialCharge::NEUTRAL:
            return 0.0;
        case PartialCharge::PARTIAL:
            return ComponentHelper::GetPartialCharge(
                atom->GetResidue(),
                atom->GetSpot(),
                atom->GetStructure());
        case PartialCharge::AMBER:
            return ComponentHelper::GetPartialCharge(
                atom->GetResidue(),
                atom->GetSpot(),
                atom->GetStructure(), true);
        default:
            Logger::Log(LogLevel::Error,
                "Invalid partial charge choice: "
                + std::to_string(static_cast<int>(m_options.partial_charge_choice)));
            break;
    }
    return 0.0;
}

void MapSimulationCommand::CalculateAtomRange()
{
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

std::unique_ptr<MapObject> MapSimulationCommand::CreateMapObject()
{
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
    return {
        std::floor(m_atom_range_minimum[0] / grid_spacing[0]) * grid_spacing[0],
        std::floor(m_atom_range_minimum[1] / grid_spacing[1]) * grid_spacing[1],
        std::floor(m_atom_range_minimum[2] / grid_spacing[2]) * grid_spacing[2]
    };
}

} // namespace rhbm_gem
