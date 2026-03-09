#pragma once

#include <array>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include <rhbm_gem/utils/FilePathHelper.hpp>

namespace ChimeraXHelper
{

struct RGB { float r = 0.2f, g = 0.6f, b = 1.0f; };

inline bool WriteCMMPoints(
    const std::vector<std::array<float,3>> & point_list,
    const std::string & path,
    float radius = 1.0f,
    RGB color = {},
    const std::string & marker_set_name = "points")
{
    std::ofstream file(path, std::ios::binary);
    if (!file) return false;

    file.setf(std::ios::fixed, std::ios::floatfield);
    file << std::setprecision(6);

    file << "<marker_set name=\"" << marker_set_name << "\">\n";
    for (size_t i = 0; i < point_list.size(); i++)
    {
        const auto & point{ point_list[i] };
        file << "  <marker id=\"" << (i + 1)
             << "\" x=\"" << point[0]
             << "\" y=\"" << point[1]
             << "\" z=\"" << point[2]
             << "\" r=\"" << color.r
             << "\" g=\"" << color.g
             << "\" b=\"" << color.b
             << "\" radius=\"" << radius
             << "\"/>\n";
    }
    file << "</marker_set>\n";
    return true;
}

inline bool WriteCMMPoints(
    const std::vector<std::array<float,3>> & point_list,
    const std::string & path,
    const std::vector<float> * radius_list,
    const std::vector<RGB> * color_list,
    float default_radius = 1.0f,
    RGB default_color = {},
    const std::string & marker_set_name = "points")
{
    std::ofstream file(path, std::ios::binary);
    if (!file) return false;

    file.setf(std::ios::fixed, std::ios::floatfield);
    file << std::setprecision(6);

    file << "<marker_set name=\"" << marker_set_name << "\">\n";
    for (size_t i = 0; i < point_list.size(); i++)
    {
        const auto & point{ point_list[i] };
        const float radius{
            (radius_list  && i < radius_list->size()) ? (*radius_list)[i] : default_radius
        };
        const RGB color{
            (color_list && i < color_list->size()) ? (*color_list)[i] : default_color
        };
        file << "  <marker id=\"" << (i + 1)
             << "\" x=\"" << point[0]
             << "\" y=\"" << point[1]
             << "\" z=\"" << point[2]
             << "\" r=\"" << color.r
             << "\" g=\"" << color.g
             << "\" b=\"" << color.b
             << "\" radius=\"" << radius
             << "\"/>\n";
    }
    file << "</marker_set>\n";
    return true;
}

} // namespace ChimeraXHelper
