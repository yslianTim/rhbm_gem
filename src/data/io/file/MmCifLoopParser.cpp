#include "internal/file/MmCifLoopParser.hpp"

#include "internal/file/MmCifTokenizer.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace rhbm_gem {

MmCifParsedDocument ParseMmCifDocumentLines(
    const std::vector<std::string> & line_list,
    const std::string & source_name)
{
    (void)source_name;

    MmCifParsedDocument parsed_document;

    size_t line_idx{ 0 };
    while (line_idx < line_list.size())
    {
        const auto & raw_line{ line_list[line_idx] };
        const auto trimmed_line{ TrimLeft(raw_line) };
        if (trimmed_line.empty() || trimmed_line.front() == '#')
        {
            ++line_idx;
            continue;
        }

        if (trimmed_line == "loop_")
        {
            ++line_idx;
            std::vector<std::string> column_name_list;
            while (line_idx < line_list.size())
            {
                const auto & header_raw_line{ line_list[line_idx] };
                const auto header_trimmed_line{ TrimLeft(header_raw_line) };
                if (header_trimmed_line.empty())
                {
                    ++line_idx;
                    continue;
                }
                if (header_trimmed_line.front() != '_') break;
                auto header_token_list{ SplitMmCifTokens(std::string{header_trimmed_line}) };
                if (!header_token_list.empty())
                {
                    column_name_list.emplace_back(header_token_list.front());
                }
                ++line_idx;
            }
            if (column_name_list.empty()) continue;

            MmCifParsedLoopCategory loop_category;
            loop_category.column_name_list = column_name_list;
            const auto expected_column_size{ loop_category.column_name_list.size() };
            std::vector<std::string> row_token_list;
            size_t row_start_line_number{ line_idx + 1 };

            while (line_idx < line_list.size())
            {
                const auto & row_raw_line{ line_list[line_idx] };
                const auto row_trimmed_line{ TrimLeft(row_raw_line) };
                if (row_trimmed_line.empty())
                {
                    ++line_idx;
                    continue;
                }
                if (row_trimmed_line.front() == '#')
                {
                    ++line_idx;
                    break;
                }
                if (row_token_list.empty() &&
                    (row_trimmed_line == "loop_" ||
                     StartsWithToken(row_trimmed_line, "data_") ||
                     row_trimmed_line.front() == '_'))
                {
                    break;
                }
                if (row_token_list.empty()) row_start_line_number = line_idx + 1;

                if (!row_raw_line.empty() && row_raw_line.front() == ';')
                {
                    std::string multiline_value{ row_raw_line.substr(1) };
                    ++line_idx;
                    auto terminated{ false };
                    while (line_idx < line_list.size())
                    {
                        const auto & multiline_raw_line{ line_list[line_idx] };
                        if (!multiline_raw_line.empty() && multiline_raw_line.front() == ';')
                        {
                            terminated = true;
                            break;
                        }
                        if (!multiline_value.empty()) multiline_value += "\n";
                        multiline_value += multiline_raw_line;
                        ++line_idx;
                    }
                    if (!terminated)
                    {
                        Logger::Log(LogLevel::Warning,
                            "ParseMmCifDocument() unterminated multiline token in loop category "
                            + BuildLoopCategoryPrefix(loop_category.column_name_list.front())
                            + " near file line " + std::to_string(row_start_line_number) + ".");
                        break;
                    }
                    row_token_list.emplace_back(std::move(multiline_value));
                    ++line_idx;
                }
                else
                {
                    auto token_list{ SplitMmCifTokens(std::string{row_trimmed_line}) };
                    row_token_list.insert(row_token_list.end(), token_list.begin(), token_list.end());
                    ++line_idx;
                }

                while (row_token_list.size() >= expected_column_size)
                {
                    std::vector<std::string> row_value_list(
                        row_token_list.begin(),
                        row_token_list.begin() + static_cast<std::ptrdiff_t>(expected_column_size));
                    loop_category.row_list.emplace_back(
                        MmCifParsedLoopRow{ std::move(row_value_list), row_start_line_number });
                    row_token_list.erase(
                        row_token_list.begin(),
                        row_token_list.begin() + static_cast<std::ptrdiff_t>(expected_column_size));
                    row_start_line_number = line_idx + 1;
                }
            }

            if (!row_token_list.empty())
            {
                Logger::Log(LogLevel::Warning,
                    "ParseMmCifDocument() loop category "
                    + BuildLoopCategoryPrefix(loop_category.column_name_list.front())
                    + " ends with incomplete row (token count = "
                    + std::to_string(row_token_list.size())
                    + ", expected = " + std::to_string(expected_column_size) + ").");
            }

            auto category_prefix{ BuildLoopCategoryPrefix(loop_category.column_name_list.front()) };
            parsed_document.loop_category_map[category_prefix].emplace_back(std::move(loop_category));
            continue;
        }

        if (trimmed_line.front() == '_')
        {
            auto token_list{ SplitMmCifTokens(std::string{trimmed_line}) };
            if (token_list.empty())
            {
                ++line_idx;
                continue;
            }

            std::string key{ token_list.front() };
            std::string value;
            if (token_list.size() >= 2)
            {
                value = token_list[1];
                parsed_document.data_item_map[key].emplace_back(std::move(value));
                ++line_idx;
                continue;
            }

            ++line_idx;
            while (line_idx < line_list.size())
            {
                const auto & value_raw_line{ line_list[line_idx] };
                const auto value_trimmed_line{ TrimLeft(value_raw_line) };
                if (value_trimmed_line.empty())
                {
                    ++line_idx;
                    continue;
                }
                if (value_trimmed_line.front() == '#') break;
                if (!value_raw_line.empty() && value_raw_line.front() == ';')
                {
                    value = value_raw_line.substr(1);
                    ++line_idx;
                    while (line_idx < line_list.size())
                    {
                        const auto & multiline_raw_line{ line_list[line_idx] };
                        if (!multiline_raw_line.empty() && multiline_raw_line.front() == ';') break;
                        if (!value.empty()) value += "\n";
                        value += multiline_raw_line;
                        ++line_idx;
                    }
                    break;
                }
                auto value_token_list{ SplitMmCifTokens(std::string{value_trimmed_line}) };
                value = value_token_list.empty() ? std::string{} : value_token_list.front();
                break;
            }
            parsed_document.data_item_map[key].emplace_back(std::move(value));
            ++line_idx;
            continue;
        }

        ++line_idx;
    }

    return parsed_document;
}

} // namespace rhbm_gem
