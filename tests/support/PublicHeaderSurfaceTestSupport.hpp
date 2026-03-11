#pragma once

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "CommandTestHelpers.hpp"

namespace contract_test_support {

inline std::vector<std::string> CollectPublicHeadersForDomain(const std::string & domain)
{
    const auto include_root{ command_test::ProjectRootPath() / "include" / "rhbm_gem" };
    const auto domain_root{ include_root / domain };
    std::vector<std::string> headers;

    for (const auto & entry : std::filesystem::recursive_directory_iterator(domain_root))
    {
        if (!entry.is_regular_file())
        {
            continue;
        }
        headers.push_back(std::filesystem::relative(entry.path(), include_root).generic_string());
    }
    std::sort(headers.begin(), headers.end());
    return headers;
}

} // namespace contract_test_support
