#pragma once
#include <filesystem>
#include <rhbm_gem/data/object/JointAnalysisResult.hpp>

namespace rhbm_gem {
// Writes saved outcomes without needing the source map or rerunning assessment.
// Throws on invalid records or output failure; does not create directories.
void WriteJointAnalysisResult(const JointAnalysisResult & result,
    const std::filesystem::path & json_path, const std::filesystem::path & csv_path);
}
