#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem {

std::vector<std::string> SplitMmCifTokens(const std::string & line);
std::string_view TrimLeft(std::string_view value);
bool StartsWithToken(std::string_view line, std::string_view token);
std::string BuildLoopCategoryPrefix(const std::string & column_name);

} // namespace rhbm_gem
