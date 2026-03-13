#include "MapSimulationCommand.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "internal/CommandDataLoader.hpp"
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/io/FileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include "workflow/DataObjectWorkflowOps.hpp"
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/math/ElectricPotential.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/core/command/OptionEnumTraits.hpp>

#include <algorithm>
#include <limits>

namespace {
constexpr std::string_view kModelKey{ "model" };
constexpr std::string_view kModelOption{ "--model" };
constexpr std::string_view kPotentialModelOption{ "--potential-model" };
constexpr std::string_view kChargeOption{ "--charge" };
constexpr std::string_view kCutoffOption{ "--cut-off" };
constexpr std::string_view kGridSpacingOption{ "--grid-spacing" };
constexpr std::string_view kBlurringWidthOption{ "--blurring-width" };
}

namespace rhbm_gem {

MapSimulationCommand::MapSimulationCommand() :
    CommandWithProfileOptions<MapSimulationCommandOptions, CommandId::MapSimulation>{},
    m_selected_atom_list{}, m_atom_charge_map{}, m_model_object{ nullptr },
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

void MapSimulationCommand::ApplyRequest(const MapSimulationRequest & request)
{
    SetThreadSize(request.common.thread_size);
    SetVerboseLevel(request.common.verbose_level);
    SetFolderPath(request.common.folder_path);
    SetModelFilePath(request.model_file_path);
    SetMapFileName(request.map_file_name);
    SetPotentialModelChoice(request.potential_model_choice);
    SetPartialChargeChoice(request.partial_charge_choice);
    SetCutoffDistance(request.cutoff_distance);
    SetGridSpacing(request.grid_spacing);
    SetBlurringWidthList(request.blurring_width_list);
}

bool MapSimulationCommand::ExecuteImpl()
{
    if (BuildDataObject() == false) return false;
    RunMapSimulation();
    return true;
}

void MapSimulationCommand::ValidateOptions()
{
    ResetPrepareIssues(kBlurringWidthOption);
    if (m_options.blurring_width_list.empty())
    {
        AddValidationError(kBlurringWidthOption,
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
    SetValidatedEnumOption(
        m_options.potential_model_choice,
        value,
        kPotentialModelOption,
        PotentialModel::FIVE_GAUS_CHARGE,
        "Potential model");
}

void MapSimulationCommand::SetPartialChargeChoice(PartialCharge value)
{
    SetValidatedEnumOption(
        m_options.partial_charge_choice,
        value,
        kChargeOption,
        PartialCharge::PARTIAL,
        "Partial charge choice");
}

void MapSimulationCommand::SetCutoffDistance(double value)
{
    SetNormalizedScalarOption(
        m_options.cutoff_distance,
        value,
        kCutoffOption,
        [](double candidate) { return candidate > 0.0; },
        5.0,
        "Cutoff distance must be positive, reset to default 5.0");
}

void MapSimulationCommand::SetModelFilePath(const std::filesystem::path & value)
{
    SetRequiredExistingPathOption(m_options.model_file_path, value, kModelOption, "Model file");
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
        kGridSpacingOption,
        [](double candidate) { return candidate > 0.0; },
        0.5,
        "Grid spacing must be positive, reset to default 0.5");
}

void MapSimulationCommand::SetBlurringWidthList(const std::string & value)
{
    MutateOptions([&]()
    {
        ResetParseIssues(kBlurringWidthOption);
        const auto parsed_list{ StringHelper::ParseListOption<double>(value, ',') };
        m_options.blurring_width_list.clear();
        m_options.blurring_width_list.reserve(parsed_list.size());
        for (const auto width : parsed_list)
        {
            if (width <= 0.0)
            {
                AddNormalizationWarning(kBlurringWidthOption,
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
        m_model_object = command_data_loader::ProcessModelFile(
            m_data_manager, m_options.model_file_path, kModelKey, "model file");
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
        WriteMap(output_file_name, *map_object);
    }
}

void MapSimulationCommand::BuildAtomList(ModelObject * model_object)
{
    if (model_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "MapSimulationCommand::BuildAtomList(): model object is null.");
        return;
    }

    SimulationAtomPreparationOptions options;
    options.partial_charge_choice = m_options.partial_charge_choice;
    options.include_unknown_atoms = true;
    auto result{ PrepareSimulationAtoms(*model_object, options) };

    m_selected_atom_list = result.atom_list;
    m_atom_charge_map = result.atom_charge_map;
    if (result.has_atom)
    {
        m_atom_range_minimum = result.range_minimum;
        m_atom_range_maximum = result.range_maximum;
        m_atom_range_minimum[0] -= static_cast<float>(m_options.cutoff_distance);
        m_atom_range_minimum[1] -= static_cast<float>(m_options.cutoff_distance);
        m_atom_range_minimum[2] -= static_cast<float>(m_options.cutoff_distance);
        m_atom_range_maximum[0] += static_cast<float>(m_options.cutoff_distance);
        m_atom_range_maximum[1] += static_cast<float>(m_options.cutoff_distance);
        m_atom_range_maximum[2] += static_cast<float>(m_options.cutoff_distance);
    }

    Logger::Log(LogLevel::Info,
        "Number of selected atoms to be simulated = "
        + std::to_string(m_selected_atom_list.size()) +" / "
        + std::to_string(model_object->GetNumberOfAtom()) + " atoms.");
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
