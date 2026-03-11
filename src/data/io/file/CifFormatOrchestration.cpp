#include "internal/io/file/CifFormat.hpp"

#include <rhbm_gem/data/object/AtomicModelDataBlock.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <sstream>
#include <stdexcept>

namespace rhbm_gem {

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
    if (m_category_parser == nullptr)
    {
        throw std::runtime_error("CifFormat category parser is not initialized.");
    }
    m_category_parser->ParseLoopBlock(data_block_prefix, m_loop_category_map, table_handler);
}

} // namespace rhbm_gem
