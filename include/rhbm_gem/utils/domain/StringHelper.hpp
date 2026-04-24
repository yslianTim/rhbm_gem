#pragma once

#include <cstddef>
#include <string>
#include <vector>
#include <sstream>
#include <cctype>
#include <stdexcept>
#include <type_traits>
#include <algorithm>

namespace rhbm_gem {
namespace string_helper {

inline void ToUpperCase(std::string & text)
{
    std::transform(
        text.begin(), text.end(), text.begin(),
        [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); }
    );
}

inline void ToLowerCase(std::string & text)
{
    std::transform(
        text.begin(), text.end(), text.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); }
    );
}

inline std::string ExtractCharAsString(const std::string & input, size_t index)
{
    if (index < input.size())
    {
        return std::string(1, input[index]);
    }
    return " ";
}

inline std::string ConvertCharArrayToString(const char * input)
{
    if (input == nullptr)
    {
        return "";
    }
    std::string result;
    while (*input != '\0')
    {
        if (*input != ' ') result += *input;
        input++;
    }
    return result;
}

inline void EraseCharFromString(std::string & input, char ch)
{
    input.erase(std::remove(input.begin(), input.end(), ch), input.end());
}

inline std::string PadWithSpaces(const std::string & text, size_t length)
{
    if (text.size() >= length)
    {
        return text;
    }
    std::string result{ text };
    result.resize(length, ' ');
    return result;
}

inline void StripCarriageReturn(std::string & line)
{
    if (!line.empty() && line.back() == '\r')
    {
        line.pop_back();
    }
}

inline std::vector<std::string> SplitStringLineFromDelimiter(
    const std::string & text, char delimiter=',')
{
    std::vector<std::string> token_list;
    std::stringstream ss(text);
    std::string token;
    while (std::getline(ss, token, delimiter))
    {
        if (!token.empty())
        {
            token_list.emplace_back(token);
        }
    }
    return token_list;
}

template <typename Type>
std::vector<Type> ParseListOption(const std::string & text, char delimiter=',')
{
    std::vector<Type> values;
    for (const auto & token : SplitStringLineFromDelimiter(text, delimiter))
    {
        if constexpr (std::is_same_v<Type, double>)
        {
            size_t pos{};
            const double value{ std::stod(token, &pos) };
            while (pos < token.size() && std::isspace(static_cast<unsigned char>(token[pos])))
            {
                ++pos;
            }
            if (pos != token.size())
            {
                throw std::invalid_argument("Invalid double: " + token);
            }
            values.emplace_back(value);
        }
        else
        {
            values.emplace_back(token);
        }
    }
    return values;
}

inline std::vector<std::string> SplitStringLineAsTokens(
    const std::string & line, size_t token_count, char delimiter='\'')
{
    std::vector<std::string> token_list;
    token_list.reserve(token_count);
    for (size_t pos = 0; pos < line.size();)
    {
        while (pos < line.size() && std::isspace(static_cast<unsigned char>(line[pos]))) ++pos;
        if (pos >= line.size()) break;
        if (line[pos] == delimiter)
        {
            ++pos;
            auto start{ pos };
            while (pos < line.size() && line[pos] != delimiter) ++pos;
            token_list.emplace_back(line.substr(start, pos - start));
            if (pos < line.size()) ++pos;
        }
        else
        {
            auto start{ pos };
            while (pos < line.size() && !std::isspace(static_cast<unsigned char>(line[pos]))) ++pos;
            token_list.emplace_back(line.substr(start, pos - start));
        }
    }
    return token_list;
}

template <typename Type>
std::string ToStringWithPrecision(const Type a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

} // namespace string_helper
} // namespace rhbm_gem
