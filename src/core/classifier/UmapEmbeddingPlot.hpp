#pragma once

#include "core/classifier/UmapEmbeddingTypes.hpp"

#include <filesystem>
#include <string>

namespace rhbm_gem::core::classifier {

struct PreparedUmapInput;

bool WriteEmbeddingPlots(
    const std::filesystem::path & output_path,
    const PreparedUmapInput & prepared,
    const UmapEmbeddingResult & first,
    const UmapEmbeddingResult & second,
    std::string & error_message);

} // namespace rhbm_gem::core::classifier
