#include "detail/CommandBase.hpp"
#include "detail/MapSampling.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/math/GridSampler.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/LocalPainter.hpp>

#include <array>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem {

class MapVisualizationCommand final : public CommandBase<MapVisualizationRequest>
{
public:
    MapVisualizationCommand();

private:
    void NormalizeAndValidateRequest(MapVisualizationRequest & request) override;
    bool ExecuteImpl(const MapVisualizationRequest & request) override;
};

namespace {

struct MapVisualizationInputs
{
    std::unique_ptr<rhbm_gem::ModelObject> model_object;
    std::unique_ptr<rhbm_gem::MapObject> map_object;
};

struct ModelAtomBondContext
{
    std::unordered_map<int, rhbm_gem::AtomObject *> atom_map;
    std::unordered_map<int, std::vector<rhbm_gem::BondObject *>> bond_map;
};

ModelAtomBondContext BuildModelAtomBondContext(rhbm_gem::ModelObject & model_object)
{
    ModelAtomBondContext result;
    const auto & selected_atom_list{ model_object.GetSelectedAtoms() };
    const auto & selected_bond_list{ model_object.GetSelectedBonds() };
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

std::filesystem::path BuildMapVisualizationOutputPath(
    const std::filesystem::path & output_dir,
    std::string_view stem,
    std::string_view extension)
{
    std::string normalized_extension{ extension };
    if (!normalized_extension.empty() && normalized_extension.front() != '.')
    {
        normalized_extension.insert(normalized_extension.begin(), '.');
    }
    return output_dir / (std::string(stem) + normalized_extension);
}

std::optional<MapVisualizationInputs> LoadMapVisualizationInputs(
    const rhbm_gem::MapVisualizationRequest & request)
{
    ScopeTimer timer("MapVisualizationCommand::BuildDataObject");
    try
    {
        MapVisualizationInputs inputs;
        inputs.model_object = rhbm_gem::ReadModel(request.model_file_path);
        inputs.model_object->SetKeyTag("model");
        inputs.map_object = rhbm_gem::ReadMap(request.map_file_path);
        inputs.map_object->SetKeyTag("map");
        return std::move(inputs);
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "MapVisualizationCommand::BuildDataObject : " + std::string(e.what()));
        return std::nullopt;
    }
}

void RunModelObjectPreprocessing(rhbm_gem::ModelObject & model_object)
{
    ScopeTimer timer("MapVisualizationCommand::RunModelObjectPreprocessing");
    model_object.SelectAllAtoms();
    model_object.SelectAllBonds();
    Logger::Log(LogLevel::Info,
        "Number of selected atom = "+ std::to_string(model_object.GetSelectedAtomCount()));
    Logger::Log(LogLevel::Info,
        "Number of selected bond = "+ std::to_string(model_object.GetSelectedBondCount()));
}

std::filesystem::path BuildOutputFilePath(
    const rhbm_gem::MapVisualizationRequest & request,
    const rhbm_gem::ModelObject & model_object,
    const std::filesystem::path & output_dir)
{
    std::string model_name;
    if (!model_object.GetPdbID().empty())
    {
        model_name = model_object.GetPdbID();
    }
    else
    {
        model_name = rhbm_gem::path_helper::GetFileName(request.model_file_path, false);
        if (model_name.empty())
        {
            model_name = "model";
        }
    }

    const std::string file_name{
        "map_slice_" + model_name + "_atom" + std::to_string(request.atom_serial_id)
    };
    return BuildMapVisualizationOutputPath(output_dir, file_name, ".pdf");
}

bool RunAtomMapValueSampling(
    const rhbm_gem::MapVisualizationRequest & request,
    rhbm_gem::MapObject & map_object,
    rhbm_gem::ModelObject & model_object,
    const std::filesystem::path & output_dir)
{
    ScopeTimer timer("MapVisualizationCommand::RunAtomMapValueSampling");
    GridSampler sampler;
    sampler.SetGridResolution(static_cast<unsigned int>(request.sampling_size));
    sampler.SetWindowSize(static_cast<float>(request.window_size));
    sampler.Print();

    auto context{ BuildModelAtomBondContext(model_object) };
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
    sampler.SetReferenceUVector({u_vector(0), u_vector(1), u_vector(2)});

    auto sampling_data_list{
        rhbm_gem::SampleMapValues(
            map_object, sampler, target_atom_position, {n_vector(0), n_vector(1), n_vector(2)})
    };
    if (sampling_data_list.empty())
    {
        Logger::Log(LogLevel::Error, "Map sampling produced no data points.");
        return false;
    }
    Eigen::MatrixXd map_value_matrix{
        Eigen::MatrixXd::Zero(request.sampling_size, request.sampling_size)
    };
    for (std::size_t index = 0; index < sampling_data_list.size(); ++index)
    {
        const std::size_t row{ index / static_cast<std::size_t>(request.sampling_size) };
        const std::size_t col{ index % static_cast<std::size_t>(request.sampling_size) };
        if (row >= static_cast<std::size_t>(request.sampling_size))
        {
            Logger::Log(LogLevel::Error,
                "Map sampling produced more values than the configured heatmap grid can hold.");
            return false;
        }
        map_value_matrix(
            static_cast<Eigen::Index>(row),
            static_cast<Eigen::Index>(col)) = static_cast<double>(sampling_data_list[index].response);
    }

    auto x_min{ target_atom_position.at(0) - 0.5 * request.window_size };
    auto x_max{ target_atom_position.at(0) + 0.5 * request.window_size };
    auto y_min{ target_atom_position.at(1) - 0.5 * request.window_size };
    auto y_max{ target_atom_position.at(1) + 0.5 * request.window_size };

    rhbm_gem::HeatmapRequest heatmap_request;
    heatmap_request.output_path = BuildOutputFilePath(request, model_object, output_dir);
    heatmap_request.values = std::move(map_value_matrix);
    heatmap_request.x_axis.title = "#font[2]{u} Position #[]{#AA}";
    heatmap_request.x_axis.range = std::make_pair(x_min, x_max);
    heatmap_request.y_axis.title = "#font[2]{v} Position #[]{#AA}";
    heatmap_request.y_axis.range = std::make_pair(y_min, y_max);
    heatmap_request.z_axis.title = "Normalized Map Value";

    const auto plot_result{ rhbm_gem::local_painter::SaveHeatmap(heatmap_request) };
    if (!plot_result.Succeeded())
    {
        Logger::Log(LogLevel::Error,
            "MapVisualizationCommand failed to render the sampled heatmap: " + plot_result.message);
        return false;
    }

    return true;
}

} // namespace

MapVisualizationCommand::MapVisualizationCommand() : CommandBase<MapVisualizationRequest>{}
{
}

void MapVisualizationCommand::NormalizeAndValidateRequest(MapVisualizationRequest & request)
{
    ValidateRequiredPath(request.model_file_path, "--model", "Model file");
    ValidateRequiredPath(request.map_file_path, "--map", "Map file");
    CoercePositiveScalar(request.atom_serial_id, "--atom-id", 1, LogLevel::Error, "Atom serial ID");
    CoercePositiveScalar(request.sampling_size, "--sampling", 100, LogLevel::Warning, "Sampling size");
    CoerceFinitePositiveScalar(request.window_size, "--window-size", 5.0, LogLevel::Error, "Window size");
}

bool MapVisualizationCommand::ExecuteImpl(const MapVisualizationRequest & request)
{
    auto inputs{ LoadMapVisualizationInputs(request) };
    if (!inputs.has_value()) return false;
    inputs->map_object->MapValueArrayNormalization();
    RunModelObjectPreprocessing(*inputs->model_object);
    return RunAtomMapValueSampling(request, *inputs->map_object, *inputs->model_object, request.output_dir);
}

namespace command_internal {

CommandResult ExecuteMapVisualizationCommand(const MapVisualizationRequest & request)
{
    MapVisualizationCommand command;
    return command.ExecuteRequest(request);
}

} // namespace command_internal

} // namespace rhbm_gem
