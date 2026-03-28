#include "internal/command/MapVisualizationCommand.hpp"
#include "internal/command/CommandCliSupport.hpp"
#include "internal/command/CommandDataSupport.hpp"
#include "internal/command/MapSampling.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/math/GridSampler.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/LocalPainter.hpp>

#include <unordered_set>
#include <tuple>
#include <vector>
#include <atomic>
#include <utility>
#include <iterator>
#include <stdexcept>
#include <limits>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace {
constexpr std::string_view kModelKey{ "model" };
constexpr std::string_view kMapKey{ "map" };
constexpr std::string_view kModelOption{ "--model" };
constexpr std::string_view kMapOption{ "--map" };
constexpr std::string_view kAtomIdOption{ "--atom-id" };
constexpr std::string_view kSamplingOption{ "--sampling" };
constexpr std::string_view kWindowSizeOption{ "--window-size" };
}

namespace rhbm_gem {

MapVisualizationCommand::MapVisualizationCommand(CommonOptionProfile profile) :
    CommandBase{ profile },
    m_model_key_tag{ kModelKey }, m_map_key_tag{ kMapKey },
    m_map_object{ nullptr }, m_model_object{ nullptr }
{
}

void BindMapVisualizationRequestOptions(
    CLI::App * command,
    MapVisualizationRequest & request)
{
    command_cli::AddPathOption(
        command,
        "-a,--model",
        [&](const std::filesystem::path & value) { request.model_file_path = value; },
        "Model file path",
        std::nullopt,
        true);
    command_cli::AddPathOption(
        command,
        "-m,--map",
        [&](const std::filesystem::path & value) { request.map_file_path = value; },
        "Map file path",
        std::nullopt,
        true);
    command_cli::AddScalarOption<int>(
        command,
        "-i,--atom-id",
        [&](int value) { request.atom_serial_id = value; },
        "Atom serial ID for visualization",
        request.atom_serial_id);
    command_cli::AddScalarOption<int>(
        command,
        "-s,--sampling",
        [&](int value) { request.sampling_size = value; },
        "Number of sampling points per atom",
        request.sampling_size);
    command_cli::AddScalarOption<double>(
        command,
        "--window-size",
        [&](double value) { request.window_size = value; },
        "Window size for sampling",
        request.window_size);
}

void MapVisualizationCommand::ApplyRequest(const MapVisualizationRequest & request)
{
    ApplyCommonRequest(request.common);
    SetModelFilePath(request.model_file_path);
    SetMapFilePath(request.map_file_path);
    SetAtomSerialID(request.atom_serial_id);
    SetSamplingSize(request.sampling_size);
    SetWindowSize(request.window_size);
}

bool MapVisualizationCommand::ExecuteImpl()
{
    if (BuildDataObject() == false) return false;
    RunMapObjectPreprocessing();
    RunModelObjectPreprocessing();
    return RunAtomMapValueSampling();
}

void MapVisualizationCommand::ResetRuntimeState()
{
    m_map_object.reset();
    m_model_object.reset();
}

void MapVisualizationCommand::SetModelFilePath(const std::filesystem::path & path)
{
    SetRequiredExistingPathOption(m_options.model_file_path, path, kModelOption, "Model file");
}

void MapVisualizationCommand::SetMapFilePath(const std::filesystem::path & path)
{
    SetRequiredExistingPathOption(m_options.map_file_path, path, kMapOption, "Map file");
}

void MapVisualizationCommand::SetAtomSerialID(int value)
{
    SetPositiveScalarOption(
        m_options.atom_serial_id,
        value,
        kAtomIdOption,
        1,
        "Atom serial ID must be positive.");
}

void MapVisualizationCommand::SetSamplingSize(int value)
{
    SetNormalizedScalarOption(
        m_options.sampling_size,
        value,
        kSamplingOption,
        [](int candidate) { return candidate > 0; },
        100,
        "Sampling size must be positive, reset to default value = 100");
}

void MapVisualizationCommand::SetWindowSize(double value)
{
    SetFinitePositiveScalarOption(
        m_options.window_size,
        value,
        kWindowSizeOption,
        5.0,
        "Window size must be a finite positive value.");
}

bool MapVisualizationCommand::BuildDataObject()
{
    ScopeTimer timer("MapVisualizationCommand::BuildDataObject");
    try
    {
        m_model_object = command_data_loader::ProcessModelFile(
            m_data_manager, m_options.model_file_path, m_model_key_tag, "model file");
        m_map_object = command_data_loader::ProcessMapFile(
            m_data_manager, m_options.map_file_path, m_map_key_tag, "map file");
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "MapVisualizationCommand::Execute() : " + std::string(e.what()));
        return false;
    }
    return true;
}

void MapVisualizationCommand::RunMapObjectPreprocessing()
{
    ScopeTimer timer("MapVisualizationCommand::RunMapObjectPreprocessing");
    NormalizeMapObject(*m_map_object);
}

void MapVisualizationCommand::RunModelObjectPreprocessing()
{
    ScopeTimer timer("MapVisualizationCommand::RunModelObjectPreprocessing");
    ModelPreparationOptions options;
    options.select_all_atoms = true;
    options.select_all_bonds = true;
    options.update_model = true;
    PrepareModelObject(*m_model_object, options);
    Logger::Log(LogLevel::Info,
        "Number of selected atom = " + std::to_string(m_model_object->GetNumberOfSelectedAtom()));
    Logger::Log(LogLevel::Info,
        "Number of selected bond = " + std::to_string(m_model_object->GetNumberOfSelectedBond()));
}

bool MapVisualizationCommand::RunAtomMapValueSampling()
{
    ScopeTimer timer("MapVisualizationCommand::RunAtomMapValueSampling");
    if (m_map_object == nullptr || m_model_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Map or model object is unavailable for visualization.");
        return false;
    }
    auto sampler{ std::make_unique<GridSampler>() };
    sampler->SetSamplingSize(static_cast<unsigned int>(m_options.sampling_size));
    sampler->SetWindowSize(static_cast<float>(m_options.window_size));
    sampler->Print();
    
    auto context{ BuildModelAtomBondContext(*m_model_object) };
    const auto & atom_map{ context.atom_map };
    const auto & bond_map{ context.bond_map };

    const auto atom_iter{ atom_map.find(m_options.atom_serial_id) };
    if (atom_iter == atom_map.end())
    {
        Logger::Log(LogLevel::Error,
            "Cannot find atom serial ID " + std::to_string(m_options.atom_serial_id)
            + " in the selected atom list.");
        return false;
    }

    const auto bond_iter{ bond_map.find(m_options.atom_serial_id) };
    if (bond_iter == bond_map.end() || bond_iter->second.empty())
    {
        Logger::Log(LogLevel::Error,
            "Cannot derive a bond context for atom serial ID "
            + std::to_string(m_options.atom_serial_id) + ".");
        return false;
    }

    auto target_atom_position{ atom_iter->second->GetPosition() };
    std::array<float, 3> reference_u_vector{ 0.0f, 0.0f, 0.0f };
    bool found_reference_bond{ false };
    for (const auto * bond : bond_iter->second)
    {
        const auto bond_vector{ bond->GetBondVector() };
        const Eigen::Map<const Eigen::Vector3f> candidate_vector(bond_vector.data());
        if (candidate_vector.norm() == 0.0f) continue;
        reference_u_vector = bond_vector;
        found_reference_bond = true;
        break;
    }
    if (!found_reference_bond)
    {
        Logger::Log(LogLevel::Error,
            "Cannot visualize atom serial ID " + std::to_string(m_options.atom_serial_id)
            + " because all available bond vectors are degenerate.");
        return false;
    }

    const Eigen::Map<const Eigen::Vector3f> eigen_u_vector(reference_u_vector.data());
    Eigen::Vector3f u_vector{ eigen_u_vector.normalized() };
    Eigen::Vector3f v_vector{ Eigen::Vector3f{0.0, 0.0, 1.0} };
    Eigen::Vector3f n_vector{ u_vector.cross(v_vector) };
    if (n_vector.norm() == 0.0f)
    {
        v_vector = Eigen::Vector3f{1.0, 0.0, 0.0};
        n_vector = u_vector.cross(v_vector);
        if (n_vector.norm() == 0.0f)
        {
            Logger::Log(LogLevel::Error,
                "Cannot construct a stable visualization axis for atom serial ID "
                + std::to_string(m_options.atom_serial_id) + ".");
            return false;
        }
    }
    sampler->SetReferenceUVector({u_vector(0), u_vector(1), u_vector(2)});

    auto sampling_data_list{
        SampleMapValues(
            *m_map_object,
            *sampler,
            target_atom_position,
            {n_vector(0), n_vector(1), n_vector(2)})
    };
    if (sampling_data_list.empty())
    {
        Logger::Log(LogLevel::Error, "Map sampling produced no data points.");
        return false;
    }
    std::vector<double> map_value_list;
    map_value_list.reserve(sampling_data_list.size());
    for (auto & [distance, value] : sampling_data_list)
    {
        map_value_list.emplace_back(static_cast<double>(value));
    }

    auto x_min{ target_atom_position.at(0) - 0.5 * m_options.window_size };
    auto x_max{ target_atom_position.at(0) + 0.5 * m_options.window_size };
    auto y_min{ target_atom_position.at(1) - 0.5 * m_options.window_size };
    auto y_max{ target_atom_position.at(1) + 0.5 * m_options.window_size };

    LocalPainter::PaintHistogram2D(
        map_value_list,
        m_options.sampling_size, x_min, x_max,
        m_options.sampling_size, y_min, y_max,
        "#font[2]{u} Position #[]{#AA}",
        "#font[2]{v} Position #[]{#AA}",
        "Normalized Map Value",
        BuildOutputFilePath().string()
    );

    return true;
}

std::filesystem::path MapVisualizationCommand::BuildOutputFilePath() const
{
    std::string model_name;
    if (m_model_object && !m_model_object->GetPdbID().empty())
    {
        model_name = m_model_object->GetPdbID();
    }
    else
    {
        model_name = FilePathHelper::GetFileName(m_options.model_file_path, false);
        if (model_name.empty())
        {
            model_name = "model";
        }
    }

    const std::string file_name{
        "map_slice_" + model_name + "_atom" + std::to_string(m_options.atom_serial_id)
    };
    return BuildOutputPath(file_name, ".pdf");
}

} // namespace rhbm_gem
