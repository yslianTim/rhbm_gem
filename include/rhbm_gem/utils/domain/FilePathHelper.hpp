#pragma once

#include <string>
#include <filesystem>

namespace rhbm_gem {
namespace path_helper {

std::string GetExtension(const std::filesystem::path & path);
std::string GetDirectory(const std::filesystem::path & path);
std::string GetFileName(const std::filesystem::path & path, bool include_extension = true);
std::string EnsureTrailingSlash(const std::filesystem::path & path);
bool EnsureFileExists(const std::filesystem::path & path, const std::string & log_prefix);

} // namespace path_helper
} // namespace rhbm_gem
