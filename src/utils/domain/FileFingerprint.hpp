#pragma once
#include <filesystem>
#include <string>

namespace rhbm_gem {
std::string FileSha256(const std::filesystem::path & path);
}
