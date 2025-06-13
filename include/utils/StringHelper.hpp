#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <sstream>

class StringHelper
{
public:
    static std::string ExtractCharAsString(const std::string & input, size_t index)
    {
        if (index < input.size())
        {
            return std::string(1, input[index]);
        }
        return " ";
    }
    
    static std::string ConvertCharArrayToString(const char * input)
    {
        std::string result;
        while (*input != '\0')
        {
            if (*input != ' ') result += *input;
            input++;
        }
        return result;
    }

    static std::vector<std::string> SpliteStringLineAsTokens(
        const std::string & line, size_t token_count, char group_delimiter='\'')
    {
        std::vector<std::string> token_list;
        token_list.reserve(token_count);
        for (size_t pos = 0; pos < line.size();)
        {
            while (pos < line.size() && std::isspace(static_cast<unsigned char>(line[pos]))) ++pos;
            if (pos >= line.size()) break;
            if (line[pos] == group_delimiter)
            {
                ++pos;
                auto start{ pos };
                while (pos < line.size() && line[pos] != group_delimiter) ++pos;
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
    static std::string ToStringWithPrecision(const Type a_value, const int n = 6)
    {
        std::ostringstream out;
        out.precision(n);
        out << std::fixed << a_value;
        return std::move(out).str();
    }

};