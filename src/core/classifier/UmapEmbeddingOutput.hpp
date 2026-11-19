#pragma once

#include "core/classifier/UmapEmbeddingTypes.hpp"

#include <filesystem>
#include <string>

namespace rhbm_gem::core::classifier {

struct PreparedUmapInput;

std::filesystem::path GetEmbeddingOutputPath(
    const std::filesystem::path & output_dir,
    const std::string & model_key_tag);

bool WriteEmbeddingOutputs(
    const std::filesystem::path & output_path,
    const PreparedUmapInput & prepared,
    const UmapEmbeddingResult & first,
    const UmapEmbeddingResult & second,
    std::string & error_message);

} // namespace rhbm_gem::core::classifier
