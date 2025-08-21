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

#include "FilePathHelper.hpp"

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

// =========================
// Quantile-tuned colormap
// =========================

// Helpers
inline float Clamp01(float x) { return x < 0.f ? 0.f : (x > 1.f ? 1.f : x); }

inline ChimeraXHelper::RGB Lerp(
    const ChimeraXHelper::RGB & a, const ChimeraXHelper::RGB & b, float t)
{
    t = Clamp01(t);
    return { a.r + (b.r - a.r) * t,
             a.g + (b.g - a.g) * t,
             a.b + (b.b - a.b) * t };
}

// Sequential palette: blue -> cyan -> yellow -> red
inline ChimeraXHelper::RGB SequentialBlueCyanYellowRed(float t)
{
    const ChimeraXHelper::RGB s0{ 0.231f, 0.298f, 0.752f }; // blue-ish
    const ChimeraXHelper::RGB s1{ 0.000f, 0.600f, 0.800f }; // cyan
    const ChimeraXHelper::RGB s2{ 0.950f, 0.900f, 0.250f }; // yellow
    const ChimeraXHelper::RGB s3{ 0.800f, 0.100f, 0.100f }; // red
    t = Clamp01(t);
    if (t < 1.0f/3.0f)      return Lerp(s0, s1, t * 3.0f);
    else if (t < 2.0f/3.0f) return Lerp(s1, s2, (t - 1.0f/3.0f) * 3.0f);
    else                    return Lerp(s2, s3, (t - 2.0f/3.0f) * 3.0f);
}

// Diverging palette centered at 0: blue -> white -> red, input x in [-1, 1]
inline ChimeraXHelper::RGB DivergingBlueWhiteRed(float x)
{
    const ChimeraXHelper::RGB blue { 0.231f, 0.298f, 0.752f };
    const ChimeraXHelper::RGB white{ 1.000f, 1.000f, 1.000f };
    const ChimeraXHelper::RGB red  { 0.865f, 0.118f, 0.118f };
    if (x < -1.f) x = -1.f;
    else if (x > 1.f) x = 1.f;
    return (x < 0.f) ? Lerp(blue, white, -x) : Lerp(white, red, x);
}

// Compute linear-interpolated quantile q in [0,1] (filters non-finite values)
inline float ComputeQuantile(const std::vector<float> & data, float q)
{
    if (std::isnan(q)) q = 0.5f;
    if (q < 0.f) q = 0.f;
    else if (q > 1.f) q = 1.f;
    std::vector<float> v;
    v.reserve(data.size());
    for (float x : data) if (std::isfinite(x)) v.push_back(x);
    if (v.empty()) return std::numeric_limits<float>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const float pos = q * (static_cast<float>(v.size() - 1));
    const size_t lo = static_cast<size_t>(std::floor(pos));
    const size_t hi = static_cast<size_t>(std::ceil(pos));
    const float t = pos - static_cast<float>(lo);
    const float a = v[lo];
    const float b = v[hi];
    return a + (b - a) * t;
}

// User options for quantile-tuned colormap
enum class Colormap { Sequential, DivergingZeroCentered };

struct QuantileColormapOptions
{
    Colormap cmap = Colormap::Sequential; // Sequential or symmetric diverging around 0
    float q_low  = 0.02f;                 // lower quantile in [0,1]
    float q_high = 0.98f;                 // upper quantile in [0,1]
    float default_radius = 1.0f;
    std::string marker_set_name = "points";
};

// Main API: map values -> colors using quantile range, then write CMM
inline bool WriteCMMPointsFromValues_QuantileColormap(
    const std::vector<std::array<float,3>> & point_list,
    const std::vector<float> & value_list,
    const std::string & path,
    const QuantileColormapOptions & opt = {})
{
    if (point_list.size() != value_list.size()) return false;

    // Gather finite values for quantiles
    std::vector<float> finite_vals;
    finite_vals.reserve(value_list.size());
    for (float v : value_list)
    {
        if (std::isfinite(v)) finite_vals.push_back(v);
    }
    if (finite_vals.empty()) return false;

    float qlo = ComputeQuantile(finite_vals, opt.q_low);
    float qhi = ComputeQuantile(finite_vals, opt.q_high);
    if (!(qhi > qlo)) qhi = qlo + 1e-6f; // avoid zero range

    std::vector<RGB> colors(value_list.size());

    if (opt.cmap == Colormap::Sequential)
    {
        const float inv{ 1.0f / (qhi - qlo) };
        for (size_t i = 0; i < value_list.size(); i++)
        {
            float t = (value_list[i] - qlo) * inv;
            t = Clamp01(t);
            colors[i] = SequentialBlueCyanYellowRed(t);
        }
    }
    else
    { // DivergingZeroCentered
        float a = std::max(std::fabs(qlo), std::fabs(qhi));
        if (!(a > 0.f)) a = 1.0f;
        const float inva = 1.0f / a;
        for (size_t i = 0; i < value_list.size(); i++)
        {
            float x = value_list[i] * inva; // map to [-1,1] (clamped in function)
            colors[i] = DivergingBlueWhiteRed(x);
        }
    }

    // Reuse existing per-point color writer
    return WriteCMMPoints(
        point_list, path,
        /*radius_list*/ nullptr,
        /*color_list*/ &colors,
        /*default_radius*/ opt.default_radius,
        /*default_color*/ {},
        /*marker_set_name*/ opt.marker_set_name);
}

} // namespace ChimeraXHelper
