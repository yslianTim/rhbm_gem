#include "MapVisualizationCommand.hpp"
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
#include "CommandRegistry.hpp"
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
rhbm_gem::CommandRegistrar<rhbm_gem::MapVisualizationCommand> registrar_map_visualization{
    "map_visualization",
    "Run map visualization"};
}

namespace rhbm_gem {

MapVisualizationCommand::MapVisualizationCommand() :
    CommandBase(), m_options{}, m_model_key_tag{"model"}, m_map_key_tag{"map"},
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
    cmd->add_option_function<std::string>("-a,--model",
        [&](const std::string & value) { SetModelFilePath(value); },
        "Model file path")->required();
    cmd->add_option_function<std::string>("-m,--map",
        [&](const std::string & value) { SetMapFilePath(value); },
        "Map file path")->required();
    cmd->add_option_function<int>("-i,--atom-id",
        [&](int value) { SetAtomSerialID(value); },
        "Atom serial ID for visualization")->default_val(m_options.atom_serial_id);
    cmd->add_option_function<int>("-s,--sampling",
        [&](int value) { SetSamplingSize(value); },
        "Number of sampling points per atom")->default_val(m_options.sampling_size);
    cmd->add_option_function<double>("--window-size",
        [&](double value) { SetWindowSize(value); },
        "Window size for sampling")->default_val(m_options.window_size);
}

bool MapVisualizationCommand::Execute()
{
    if (BuildDataObject() == false) return false;
    RunMapObjectPreprocessing();
    RunModelObjectPreprocessing();
    RunAtomMapValueSampling();

    return true;
}

void MapVisualizationCommand::SetModelFilePath(const std::filesystem::path & path)
{
    m_options.model_file_path = path;
    if (!FilePathHelper::EnsureFileExists(m_options.model_file_path, "Model file"))
    {
        Logger::Log(LogLevel::Error,
            "Model file does not exist: " + m_options.model_file_path.string());
        m_valiate_options = false;
    }
}

void MapVisualizationCommand::SetMapFilePath(const std::filesystem::path & path)
{
    m_options.map_file_path = path;
    if (!FilePathHelper::EnsureFileExists(m_options.map_file_path, "Map file"))
    {
        Logger::Log(LogLevel::Error,
            "Map file does not exist: " + m_options.map_file_path.string());
        m_valiate_options = false;
    }
}

void MapVisualizationCommand::SetSamplingSize(int value)
{
    m_options.sampling_size = value;
    if (m_options.sampling_size <= 0)
    {
        Logger::Log(LogLevel::Warning,
            "Sampling size must be positive, reset to default value = 1500");
        m_options.sampling_size = 1500;
    }
}

void MapVisualizationCommand::SetWindowSize(double value)
{
    m_options.window_size = value;
}

bool MapVisualizationCommand::BuildDataObject()
{
    ScopeTimer timer("MapVisualizationCommand::BuildDataObject");
    auto data_manager{ GetDataManagerPtr() };
    data_manager->SetDatabaseManager(m_options.database_path);
    try
    {
        data_manager->ProcessFile(m_options.model_file_path, m_model_key_tag);
        data_manager->ProcessFile(m_options.map_file_path, m_map_key_tag);
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
    auto data_manager{ GetDataManagerPtr() };
    m_map_object = data_manager->GetTypedDataObject<MapObject>(m_map_key_tag);
    m_map_object->MapValueArrayNormalization();
}

void MapVisualizationCommand::RunModelObjectPreprocessing()
{
    ScopeTimer timer("MapVisualizationCommand::RunModelObjectPreprocessing");
    auto data_manager{ GetDataManagerPtr() };
    m_model_object = data_manager->GetTypedDataObject<ModelObject>(m_model_key_tag);
    for (auto & atom : m_model_object->GetAtomList()) atom->SetSelectedFlag(true);
    for (auto & bond : m_model_object->GetBondList()) bond->SetSelectedFlag(true);
    m_model_object->Update();
    Logger::Log(LogLevel::Info,
        "Number of selected atom = " + std::to_string(m_model_object->GetNumberOfSelectedAtom()));
    Logger::Log(LogLevel::Info,
        "Number of selected bond = " + std::to_string(m_model_object->GetNumberOfSelectedBond()));
}

void MapVisualizationCommand::RunAtomMapValueSampling()
{
    ScopeTimer timer("MapVisualizationCommand::RunAtomMapValueSampling");
    if (m_map_object == nullptr) return;
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
    }
    auto target_atom_position{ atom_map.at(m_options.atom_serial_id)->GetPosition() };
    auto reference_u_vector{ bond_map.at(m_options.atom_serial_id).at(1)->GetBondVector() };
    bond_map.at(m_options.atom_serial_id).at(1)->Display();
    const Eigen::Map<const Eigen::Vector3f> eigen_u_vector(reference_u_vector.data());
    Eigen::Vector3f u_vector{ eigen_u_vector.normalized() };
    Eigen::Vector3f v_vector{ Eigen::Vector3f{0.0, 0.0, 1.0} };
    Eigen::Vector3f n_vector{ u_vector.cross(v_vector) };
    sampler->SetReferenceUVector({u_vector(0), u_vector(1), u_vector(2)});

    MapInterpolationVisitor interpolation_visitor{ sampler.get() };
    
    interpolation_visitor.SetPosition(target_atom_position);
    interpolation_visitor.SetAxisVector({n_vector(0), n_vector(1), n_vector(2)});
    m_map_object->Accept(&interpolation_visitor);
    auto & sampling_data_list{ interpolation_visitor.GetSamplingDataList() };
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
        "/Users/yslian/Downloads/test_2d_slice_map.pdf"
    );

}

} // namespace rhbm_gem
