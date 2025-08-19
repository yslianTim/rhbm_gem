#pragma once

#include <array>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace chimerax {

struct RGB { float r = 0.2f, g = 0.6f, b = 1.0f; };

// ========== 基礎：寫 CMM（整體同色同半徑） ==========
inline bool write_cmm_points(const std::vector<std::array<float,3>>& pts,
                             const std::string& path,
                             float radius = 1.0f,          // Å
                             RGB color = {},               // 0..1
                             const std::string& marker_set_name = "points")
{
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;

    f.setf(std::ios::fixed, std::ios::floatfield);
    f << std::setprecision(6);

    f << "<marker_set name=\"" << marker_set_name << "\">\n";
    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& p = pts[i];
        f << "  <marker id=\"" << (i + 1)
          << "\" x=\"" << p[0]
          << "\" y=\"" << p[1]
          << "\" z=\"" << p[2]
          << "\" r=\"" << color.r
          << "\" g=\"" << color.g
          << "\" b=\"" << color.b
          << "\" radius=\"" << radius
          << "\"/>\n";
    }
    f << "</marker_set>\n";
    return true;
}

// ========== 進階：寫 CMM（每點可自訂半徑/顏色） ==========
// radii/colors 可為 nullptr；存在時會按索引套用（長度不足則落回預設）
inline bool write_cmm_points(const std::vector<std::array<float,3>>& pts,
                             const std::string& path,
                             const std::vector<float>* radii,
                             const std::vector<RGB>* colors,
                             float default_radius = 1.0f,
                             RGB   default_color  = {},
                             const std::string& marker_set_name = "points")
{
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;

    f.setf(std::ios::fixed, std::ios::floatfield);
    f << std::setprecision(6);

    f << "<marker_set name=\"" << marker_set_name << "\">\n";
    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& p = pts[i];
        const float rad = (radii  && i < radii->size())  ? (*radii)[i] : default_radius;
        const RGB   col = (colors && i < colors->size()) ? (*colors)[i] : default_color;
        f << "  <marker id=\"" << (i + 1)
          << "\" x=\"" << p[0]
          << "\" y=\"" << p[1]
          << "\" z=\"" << p[2]
          << "\" r=\"" << col.r
          << "\" g=\"" << col.g
          << "\" b=\"" << col.b
          << "\" radius=\"" << rad
          << "\"/>\n";
    }
    f << "</marker_set>\n";
    return true;
}

// ========== XYZ ==========
inline bool write_xyz_points(const std::vector<std::array<float,3>>& pts,
                             const std::string& path,
                             const std::string& element = "C", // 當作虛擬原子來畫球
                             const std::string& comment = "points")
{
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;

    f.setf(std::ios::fixed, std::ios::floatfield);
    f << std::setprecision(6);

    f << pts.size() << "\n" << comment << "\n";
    for (const auto& p : pts) {
        f << element << ' ' << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
    }
    return true;
}

// ========== 便利函式：依副檔名自動選擇 ==========
inline bool write_points_auto(const std::vector<std::array<float,3>>& pts,
                              const std::string& path,
                              float radius = 1.0f,
                              RGB color = {},
                              const std::string& marker_set_name = "points",
                              const std::string& element = "C",
                              const std::string& comment = "points")
{
    auto ends_with_ci = [](const std::string& s, const std::string& suf) {
        if (s.size() < suf.size()) return false;
        size_t off = s.size() - suf.size();
        for (size_t i = 0; i < suf.size(); ++i) {
            char a = static_cast<char>(std::tolower(static_cast<unsigned char>(s[off + i])));
            char b = static_cast<char>(std::tolower(static_cast<unsigned char>(suf[i])));
            if (a != b) return false;
        }
        return true;
    };

    if (ends_with_ci(path, ".cmm")) {
        return write_cmm_points(pts, path, radius, color, marker_set_name);
    } else if (ends_with_ci(path, ".xyz")) {
        return write_xyz_points(pts, path, element, comment);
    }
    return false; // 不支援的副檔名
}

} // namespace chimerax