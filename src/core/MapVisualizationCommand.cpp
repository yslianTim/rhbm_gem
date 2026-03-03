#include "MapVisualizationCommand.hpp"
#include "CommandOptionBinding.hpp"
#include "DataObjectManager.hpp"
#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "MapObject.hpp"
#include "ModelObject.hpp"
#include "MapInterpolationVisitor.hpp"
#include "ScopeTimer.hpp"
#include "FilePathHelper.hpp"
#include "ChemicalDataHelper.hpp"
#include "GridSampler.hpp"
#include "Logger.hpp"
#include "LocalPainter.hpp"

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
}

namespace rhbm_gem {

MapVisualizationCommand::MapVisualizationCommand() :
    CommandBase(), m_options{}, m_model_key_tag{ kModelKey }, m_map_key_tag{ kMapKey },
    m_map_object{ nullptr }, m_model_object{ nullptr }
{
}

MapVisualizationCommand::~MapVisualizationCommand()
{
    m_map_object.reset();
    m_model_object.reset();
}

void MapVisualizationCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    command_cli::AddPathOption(
        cmd, "-a,--model",
        [&](const std::filesystem::path & value) { SetModelFilePath(value); },
        "Model file path",
        std::nullopt,
        true);
    command_cli::AddPathOption(
        cmd, "-m,--map",
        [&](const std::filesystem::path & value) { SetMapFilePath(value); },
        "Map file path",
        std::nullopt,
        true);
    command_cli::AddScalarOption<int>(
        cmd, "-i,--atom-id",
        [&](int value) { SetAtomSerialID(value); },
        "Atom serial ID for visualization",
        m_options.atom_serial_id);
    command_cli::AddScalarOption<int>(
        cmd, "-s,--sampling",
        [&](int value) { SetSamplingSize(value); },
        "Number of sampling points per atom",
        m_options.sampling_size);
    command_cli::AddScalarOption<double>(
        cmd, "--window-size",
        [&](double value) { SetWindowSize(value); },
        "Window size for sampling",
        m_options.window_size);
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
    SetRequiredExistingPathOption(m_options.model_file_path, path, "--model", "Model file");
}

void MapVisualizationCommand::SetMapFilePath(const std::filesystem::path & path)
{
    SetRequiredExistingPathOption(m_options.map_file_path, path, "--map", "Map file");
}

void MapVisualizationCommand::SetAtomSerialID(int value)
{
    SetPositiveScalarOption(
        m_options.atom_serial_id,
        value,
        "--atom-id",
        1,
        "Atom serial ID must be positive.");
}

void MapVisualizationCommand::SetSamplingSize(int value)
{
    SetNormalizedScalarOption(
        m_options.sampling_size,
        value,
        "--sampling",
        [](int candidate) { return candidate > 0; },
        100,
        "Sampling size must be positive, reset to default value = 100");
}

void MapVisualizationCommand::SetWindowSize(double value)
{
    SetFinitePositiveScalarOption(
        m_options.window_size,
        value,
        "--window-size",
        5.0,
        "Window size must be a finite positive value.");
}

bool MapVisualizationCommand::BuildDataObject()
{
    ScopeTimer timer("MapVisualizationCommand::BuildDataObject");
    try
    {
        m_model_object = ProcessTypedFile<ModelObject>(
            m_options.model_file_path, m_model_key_tag, "model file");
        m_map_object = ProcessTypedFile<MapObject>(
            m_options.map_file_path, m_map_key_tag, "map file");
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
    m_map_object->MapValueArrayNormalization();
}

void MapVisualizationCommand::RunModelObjectPreprocessing()
{
    ScopeTimer timer("MapVisualizationCommand::RunModelObjectPreprocessing");
    for (auto & atom : m_model_object->GetAtomList()) atom->SetSelectedFlag(true);
    for (auto & bond : m_model_object->GetBondList()) bond->SetSelectedFlag(true);
    m_model_object->Update();
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
    
    const auto & atom_list{ m_model_object->GetSelectedAtomList() };
    const auto & bond_list{ m_model_object->GetSelectedBondList() };
    std::unordered_map<int, AtomObject*> atom_map;
    std::unordered_map<int, std::vector<BondObject*>> bond_map;
    for (auto & atom : atom_list)
    {
        atom_map.emplace(atom->GetSerialID(), atom);
    }
    for (auto & bond : bond_list)
    {
        bond_map[bond->GetAtomSerialID1()].emplace_back(bond);
        bond_map[bond->GetAtomSerialID2()].emplace_back(bond);
    }

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

    MapInterpolationVisitor interpolation_visitor{ sampler.get() };
    
    interpolation_visitor.SetPosition(target_atom_position);
    interpolation_visitor.SetAxisVector({n_vector(0), n_vector(1), n_vector(2)});
    m_map_object->Accept(&interpolation_visitor);
    auto & sampling_data_list{ interpolation_visitor.GetSamplingDataList() };
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
