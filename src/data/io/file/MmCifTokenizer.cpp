#include "internal/file/MmCifTokenizer.hpp"

#include <cctype>

namespace rhbm_gem {

std::vector<std::string> SplitMmCifTokens(const std::string & line)
{
    std::vector<std::string> token_list;
    std::string current_token;
    enum class State
    {
        IN_SPACE,
        IN_UNQUOTED,
        IN_SINGLE_QUOTE,
        IN_DOUBLE_QUOTE
    };
    State state{ State::IN_SPACE };

    auto flush_token{
        [&]()
        {
            if (current_token.empty()) return;
            token_list.emplace_back(std::move(current_token));
            current_token.clear();
        }
    };

    for (size_t pos = 0; pos < line.size(); ++pos)
    {
        const char current_char{ line[pos] };
        switch (state)
        {
        case State::IN_SPACE:
            if (std::isspace(static_cast<unsigned char>(current_char)))
            {
                continue;
            }
            if (current_char == '\'')
            {
                state = State::IN_SINGLE_QUOTE;
                continue;
            }
            if (current_char == '"')
            {
                state = State::IN_DOUBLE_QUOTE;
                continue;
            }
            current_token.push_back(current_char);
            state = State::IN_UNQUOTED;
            continue;
        case State::IN_UNQUOTED:
            if (std::isspace(static_cast<unsigned char>(current_char)))
            {
                flush_token();
                state = State::IN_SPACE;
                continue;
            }
            current_token.push_back(current_char);
            continue;
        case State::IN_SINGLE_QUOTE:
            if (current_char == '\'')
            {
                flush_token();
                state = State::IN_SPACE;
                continue;
            }
            current_token.push_back(current_char);
            continue;
        case State::IN_DOUBLE_QUOTE:
            if (current_char == '"')
            {
                flush_token();
                state = State::IN_SPACE;
                continue;
            }
            current_token.push_back(current_char);
            continue;
        }
    }
    flush_token();
    return token_list;
}

std::string_view TrimLeft(std::string_view value)
{
    size_t pos{ 0 };
    while (pos < value.size() && std::isspace(static_cast<unsigned char>(value[pos])))
    {
        ++pos;
    }
    return value.substr(pos);
}

bool StartsWithToken(std::string_view line, std::string_view token)
{
    return line.rfind(token, 0) == 0;
}

std::string BuildLoopCategoryPrefix(const std::string & column_name)
{
    auto dot_pos{ column_name.find('.') };
    if (dot_pos == std::string::npos) return column_name;
    return column_name.substr(0, dot_pos + 1);
}

} // namespace rhbm_gem
