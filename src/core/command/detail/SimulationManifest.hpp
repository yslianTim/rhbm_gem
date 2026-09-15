#pragma once

#include "MapSimulation.hpp"

#include <filesystem>
#include <string>

namespace rhbm_gem::core::simulation {

struct SimulationSource
{
    std::string pdb_id;
    std::string model_sha256;
};

std::string FileSha256(const std::filesystem::path & path);
void WriteSimulationArtifacts(const std::filesystem::path & output, const MapObject & map,
    const SimulationAtomPreparationResult & atoms, const MapSimulationRequest & request,
    double blurring_width, int actual_job_count, const SimulationSource & source);

} // namespace rhbm_gem::core::simulation
