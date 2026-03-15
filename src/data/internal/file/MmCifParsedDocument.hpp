#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace rhbm_gem {

struct MmCifParsedLoopRow
{
    std::vector<std::string> token_list;
    size_t line_number;
};

struct MmCifParsedLoopCategory
{
    std::vector<std::string> column_name_list;
    std::vector<MmCifParsedLoopRow> row_list;
};

struct MmCifParsedDocument
{
    std::unordered_map<std::string, std::vector<MmCifParsedLoopCategory>> loop_category_map;
    std::unordered_map<std::string, std::vector<std::string>> data_item_map;
};

} // namespace rhbm_gem
