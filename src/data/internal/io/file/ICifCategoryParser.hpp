#pragma once

#include "MmCifParsedDocument.hpp"

#include <functional>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace rhbm_gem {

using CifColumnIndexMap = std::map<std::string, size_t, std::less<>>;
using CifLoopRowHandler = std::function<void(
    const CifColumnIndexMap &,
    const std::vector<std::string> &)>;

class ICifCategoryParser
{
public:
    virtual ~ICifCategoryParser() = default;

    virtual void ParseLoopBlock(
        std::string_view data_block_prefix,
        const std::unordered_map<std::string, std::vector<MmCifParsedLoopCategory>> & loop_category_map,
        const CifLoopRowHandler & row_handler) const = 0;
};

class MmCifCategoryParser final : public ICifCategoryParser
{
public:
    void ParseLoopBlock(
        std::string_view data_block_prefix,
        const std::unordered_map<std::string, std::vector<MmCifParsedLoopCategory>> & loop_category_map,
        const CifLoopRowHandler & row_handler) const override;
};

} // namespace rhbm_gem
