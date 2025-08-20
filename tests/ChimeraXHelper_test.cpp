#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <array>

#include "ChimeraXHelper.hpp"

TEST(ChimeraXHelperTest, WritePointsCmm)
{
    using ChimeraXHelper::RGB;
    std::vector<std::array<float,3>> pts{{{1.0f,2.0f,3.0f},{4.0f,5.0f,6.0f}}};
    std::filesystem::path tmp = std::filesystem::temp_directory_path()/"chimera_test.cmm";
    ASSERT_TRUE(ChimeraXHelper::WriteCMMPoints(pts, tmp.string(), 2.0f, RGB{0.1f,0.2f,0.3f}, "test"));
    std::ifstream ifs(tmp);
    ASSERT_TRUE(ifs.good());
    std::string content((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    std::string expected =
        "<marker_set name=\"test\">\n"
        "  <marker id=\"1\" x=\"1.000000\" y=\"2.000000\" z=\"3.000000\" r=\"0.100000\" g=\"0.200000\" b=\"0.300000\" radius=\"2.000000\"/>\n"
        "  <marker id=\"2\" x=\"4.000000\" y=\"5.000000\" z=\"6.000000\" r=\"0.100000\" g=\"0.200000\" b=\"0.300000\" radius=\"2.000000\"/>\n"
        "</marker_set>\n";
    EXPECT_EQ(expected, content);
    std::filesystem::remove(tmp);
}
