#pragma once

#include "MmCifParsedDocument.hpp"

#include <string>
#include <vector>

namespace rhbm_gem {

MmCifParsedDocument ParseMmCifDocumentLines(
    const std::vector<std::string> & line_list,
    const std::string & source_name);

} // namespace rhbm_gem
