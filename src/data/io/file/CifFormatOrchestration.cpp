#include "internal/file/CifFormat.hpp"

#include "internal/file/AtomicModelDataBlock.hpp"
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <sstream>
#include <stdexcept>

namespace rhbm_gem {

namespace {

std::string BuildMmCifTokenPreview(
    const std::vector<std::string> & token_list,
    size_t max_items = 8)
{
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < token_list.size() && i < max_items; ++i)
    {
        if (i > 0) oss << ", ";
        oss << token_list[i];
    }
    if (token_list.size() > max_items) oss << ", ...";
    oss << "]";
    return oss.str();
}

void ParseMmCifLoopCategory(
    std::string_view data_block_prefix,
    const std::unordered_map<std::string, std::vector<MmCifParsedLoopCategory>> & loop_category_map,
    const CifLoopRowHandler & row_handler)
{
    auto category_iter{ loop_category_map.find(std::string{data_block_prefix}) };
    if (category_iter == loop_category_map.end()) return;

    for (const auto & category : category_iter->second)
    {
        CifColumnIndexMap column_index_map;
        for (size_t i = 0; i < category.column_name_list.size(); ++i)
        {
            std::string short_name{ category.column_name_list[i] };
            if (short_name.rfind(data_block_prefix, 0) == 0)
            {
                short_name = short_name.substr(data_block_prefix.size());
            }
            column_index_map[short_name] = i;
        }

        size_t loop_row_number{ 0 };
        for (const auto & row : category.row_list)
        {
            ++loop_row_number;
            if (row.token_list.size() != column_index_map.size())
            {
                Logger::Log(
                    LogLevel::Warning,
                    "ParseLoopBlock(" + std::string(data_block_prefix)
                        + ") row " + std::to_string(loop_row_number)
                        + " has token count " + std::to_string(row.token_list.size())
                        + " but expects " + std::to_string(column_index_map.size())
                        + " at file line " + std::to_string(row.line_number)
                        + ". tokens = " + BuildMmCifTokenPreview(row.token_list));
                continue;
            }
            try
            {
                row_handler(column_index_map, row.token_list);
            }
            catch (const std::exception & ex)
            {
                Logger::Log(
                    LogLevel::Warning,
                    "ParseLoopBlock(" + std::string(data_block_prefix)
                        + ") row " + std::to_string(loop_row_number)
                        + " failed at file line " + std::to_string(row.line_number)
                        + ": " + std::string(ex.what())
                        + ". tokens = " + BuildMmCifTokenPreview(row.token_list));
            }
        }
    }
}

} // namespace

void CifFormat::LogHeaderSummary() const
{
    std::ostringstream oss;
    oss << "CIF Header Information:\n";
    oss << "#Entities = " << m_data_block->GetEntityTypeMap().size() << "\n";
    for (const auto & [entity_id, chain_id] : m_data_block->GetChainIDListMap())
    {
        oss << "[" << entity_id << "] : ";
        for (size_t i = 0; i < chain_id.size(); i++)
        {
            oss << chain_id.at(i);
            if (i < chain_id.size() - 1) oss << ",";
        }
        oss << "\n";
    }

    const auto element_size{ m_data_block->GetElementTypeList().size() };
    oss << "#Elementry types = " << element_size << "\n";
    oss << "Element type list : ";
    size_t count{ 0 };
    for (const auto element : m_data_block->GetElementTypeList())
    {
        oss << ChemicalDataHelper::GetLabel(element);
        if (count < element_size - 1) oss << ",";
        count++;
    }
    oss << "\n";
    Logger::Log(LogLevel::Info, oss.str());
}

void CifFormat::ParseLoopBlock(
    std::string_view data_block_prefix,
    const CifLoopRowHandler & table_handler)
{
    ParseMmCifLoopCategory(data_block_prefix, m_loop_category_map, table_handler);
}

} // namespace rhbm_gem
