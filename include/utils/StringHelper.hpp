#pragma once

#include <iostream>
#include <string>
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

    template <typename Type>
    static std::string ToStringWithPrecision(const Type a_value, const int n = 6)
    {
        std::ostringstream out;
        out.precision(n);
        out << std::fixed << a_value;
        return std::move(out).str();
    }

};