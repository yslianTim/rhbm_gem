#pragma once

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

namespace command_test {

inline std::filesystem::path ProjectRootPath()
{
    return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
}

inline std::filesystem::path TestDataPath(const std::string & file_name)
{
    return std::filesystem::path(__FILE__).parent_path().parent_path() / "fixtures" / file_name;
}

inline std::filesystem::path MakeUniqueTempDir(const std::string & prefix)
{
    const auto timestamp{
        static_cast<std::uint64_t>(
            std::chrono::steady_clock::now().time_since_epoch().count())
    };
    const auto path{
        std::filesystem::temp_directory_path()
        / (prefix + "_" + std::to_string(timestamp))
    };
    std::filesystem::create_directories(path);
    return path;
}

class ScopedTempDir
{
public:
    explicit ScopedTempDir(const std::string & prefix) :
        m_path{ MakeUniqueTempDir(prefix) }
    {
    }

    ~ScopedTempDir()
    {
        std::error_code error_code;
        std::filesystem::remove_all(m_path, error_code);
    }

    const std::filesystem::path & path() const
    {
        return m_path;
    }

private:
    std::filesystem::path m_path;
};

inline void SeedSavedModel(
    const std::filesystem::path & database_path,
    const std::filesystem::path & model_path,
    const std::string & saved_key,
    const std::string & pdb_id)
{
    rhbm_gem::DataRepository repository{ database_path };
    auto model{ rhbm_gem::ReadModel(model_path) };
    model->SetKeyTag("model");
    model->SetPdbID(pdb_id);
    auto analysis{ model->EditAnalysis() };
    for (auto & atom : model->GetAtomList())
    {
        analysis.EnsureAtomLocalPotential(*atom);
    }
    repository.SaveModel(*model, saved_key);
}

inline std::filesystem::path GenerateMapFile(
    const std::filesystem::path & output_dir,
    const std::filesystem::path & model_path,
    const std::string & map_name = "sim_map",
    const std::string & blurring_widths = "1.0")
{
    rhbm_gem::MapSimulationRequest request{};
    request.output_dir = output_dir;
    request.model_file_path = model_path;
    request.map_file_name = map_name;
    request.blurring_width_list = rhbm_gem::string_helper::ParseListOption<double>(blurring_widths);

    const auto result{ rhbm_gem::RunMapSimulation(request) };
    if (!result.succeeded)
    {
        throw std::runtime_error("Failed to generate map fixture.");
    }

    for (const auto & entry : std::filesystem::directory_iterator(output_dir))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".map")
        {
            return entry.path();
        }
    }

    throw std::runtime_error("Generated map fixture was not found.");
}

inline std::size_t CountFilesWithExtension(
    const std::filesystem::path & directory,
    const std::string & extension)
{
    if (!std::filesystem::exists(directory)) return 0;

    std::size_t count{ 0 };
    for (const auto & entry : std::filesystem::directory_iterator(directory))
    {
        if (entry.is_regular_file() && entry.path().extension() == extension)
        {
            ++count;
        }
    }
    return count;
}

} // namespace command_test
