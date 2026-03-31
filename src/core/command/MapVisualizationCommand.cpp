#include "internal/command/MapVisualizationCommand.hpp"
#include "internal/command/MapSampling.hpp"
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
#include <unordered_map>
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

template <typename TypedDataObject>
std::shared_ptr<TypedDataObject> LoadTypedObjectFromFile(
    rhbm_gem::DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label,
    std::string_view object_type_name)
{
    try
    {
        data_manager.ProcessFile(path, std::string(key_tag));
        return data_manager.GetTypedDataObject<TypedDataObject>(std::string(key_tag));
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to process " + std::string(label) + " from '" + path.string()
            + "' as " + std::string(object_type_name) + ": " + ex.what());
    }
}

struct ModelAtomBondContext
{
    std::unordered_map<int, rhbm_gem::AtomObject *> atom_map;
    std::unordered_map<int, std::vector<rhbm_gem::BondObject *>> bond_map;
};

void SelectAllAtoms(rhbm_gem::ModelObject & model_object)
{
    for (auto & atom : model_object.GetAtomList())
    {
        atom->SetSelectedFlag(true);
    }
}

void SelectAllBonds(rhbm_gem::ModelObject & model_object)
{
    for (auto & bond : model_object.GetBondList())
    {
        bond->SetSelectedFlag(true);
    }
}

void PrepareModelForVisualization(rhbm_gem::ModelObject & model_object)
{
    SelectAllAtoms(model_object);
    SelectAllBonds(model_object);
    model_object.Update();
}

ModelAtomBondContext BuildModelAtomBondContext(rhbm_gem::ModelObject & model_object)
{
    ModelAtomBondContext result;
    const auto & selected_atom_list{ model_object.GetSelectedAtomList() };
    const auto & selected_bond_list{ model_object.GetSelectedBondList() };
    result.atom_map.reserve(selected_atom_list.size());
    result.bond_map.reserve(selected_atom_list.size());

    for (auto * atom : selected_atom_list)
    {
        result.atom_map.emplace(atom->GetSerialID(), atom);
    }
    for (auto * bond : selected_bond_list)
    {
        result.bond_map[bond->GetAtomSerialID1()].emplace_back(bond);
        result.bond_map[bond->GetAtomSerialID2()].emplace_back(bond);
    }
    return result;
}
}

namespace rhbm_gem {

MapVisualizationCommand::MapVisualizationCommand() :
    CommandWithRequest<MapVisualizationRequest>{},
    m_model_key_tag{ kModelKey }, m_map_key_tag{ kMapKey },
    m_map_object{ nullptr }, m_model_object{ nullptr }
{
}

void MapVisualizationCommand::NormalizeRequest()
{
    auto & request{ MutableRequest() };
    ValidateRequiredPath(
        request.model_file_path,
        kModelOption,
        "Model file");
    ValidateRequiredPath(
        request.map_file_path,
        kMapOption,
        "Map file");
    CoerceScalar(
        request.atom_serial_id,
        kAtomIdOption,
        [](int candidate) { return candidate > 0; },
        1,
        LogLevel::Error,
        "Atom serial ID must be positive.");
    CoerceScalar(
        request.sampling_size,
        kSamplingOption,
        [](int candidate) { return candidate > 0; },
        100,
        LogLevel::Warning,
        "Sampling size must be positive, reset to default value = 100");
    CoerceScalar(
        request.window_size,
        kWindowSizeOption,
        [](double candidate)
        {
            return std::isfinite(candidate) && candidate > 0.0;
        },
        5.0,
        LogLevel::Error,
        "Window size must be a finite positive value.");
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

bool MapVisualizationCommand::BuildDataObject()
{
    const auto & request{ RequestOptions() };
    ScopeTimer timer("MapVisualizationCommand::BuildDataObject");
    try
    {
        m_model_object = LoadTypedObjectFromFile<ModelObject>(
            m_data_manager, request.model_file_path, m_model_key_tag, "model file", "ModelObject");
        m_map_object = LoadTypedObjectFromFile<MapObject>(
            m_data_manager, request.map_file_path, m_map_key_tag, "map file", "MapObject");
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
    PrepareModelForVisualization(*m_model_object);
    Logger::Log(LogLevel::Info,
        "Number of selected atom = " + std::to_string(m_model_object->GetNumberOfSelectedAtom()));
    Logger::Log(LogLevel::Info,
        "Number of selected bond = " + std::to_string(m_model_object->GetNumberOfSelectedBond()));
}

bool MapVisualizationCommand::RunAtomMapValueSampling()
{
    const auto & request{ RequestOptions() };
    ScopeTimer timer("MapVisualizationCommand::RunAtomMapValueSampling");
    if (m_map_object == nullptr || m_model_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Map or model object is unavailable for visualization.");
        return false;
    }
    auto sampler{ std::make_unique<GridSampler>() };
    sampler->SetSamplingSize(static_cast<unsigned int>(request.sampling_size));
    sampler->SetWindowSize(static_cast<float>(request.window_size));
    sampler->Print();
    
    auto context{ BuildModelAtomBondContext(*m_model_object) };
    const auto & atom_map{ context.atom_map };
    const auto & bond_map{ context.bond_map };

    const auto atom_iter{ atom_map.find(request.atom_serial_id) };
    if (atom_iter == atom_map.end())
    {
        Logger::Log(LogLevel::Error,
            "Cannot find atom serial ID " + std::to_string(request.atom_serial_id)
            + " in the selected atom list.");
        return false;
    }

    const auto bond_iter{ bond_map.find(request.atom_serial_id) };
    if (bond_iter == bond_map.end() || bond_iter->second.empty())
    {
        Logger::Log(LogLevel::Error,
            "Cannot derive a bond context for atom serial ID "
            + std::to_string(request.atom_serial_id) + ".");
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
            "Cannot visualize atom serial ID " + std::to_string(request.atom_serial_id)
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
                + std::to_string(request.atom_serial_id) + ".");
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

    auto x_min{ target_atom_position.at(0) - 0.5 * request.window_size };
    auto x_max{ target_atom_position.at(0) + 0.5 * request.window_size };
    auto y_min{ target_atom_position.at(1) - 0.5 * request.window_size };
    auto y_max{ target_atom_position.at(1) + 0.5 * request.window_size };

    LocalPainter::PaintHistogram2D(
        map_value_list,
        request.sampling_size, x_min, x_max,
        request.sampling_size, y_min, y_max,
        "#font[2]{u} Position #[]{#AA}",
        "#font[2]{v} Position #[]{#AA}",
        "Normalized Map Value",
        BuildOutputFilePath().string()
    );

    return true;
}

std::filesystem::path MapVisualizationCommand::BuildOutputFilePath() const
{
    const auto & request{ RequestOptions() };
    std::string model_name;
    if (m_model_object && !m_model_object->GetPdbID().empty())
    {
        model_name = m_model_object->GetPdbID();
    }
    else
    {
        model_name = FilePathHelper::GetFileName(request.model_file_path, false);
        if (model_name.empty())
        {
            model_name = "model";
        }
    }

    const std::string file_name{
        "map_slice_" + model_name + "_atom" + std::to_string(request.atom_serial_id)
    };
    return BuildOutputPath(file_name, ".pdf");
}

} // namespace rhbm_gem
