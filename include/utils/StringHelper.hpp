#pragma once

#include <iostream>
#include <string>

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

};