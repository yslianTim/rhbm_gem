#include <gtest/gtest.h>
#include <filesystem>
#include <memory>
#include <array>
#include <vector>

#include <rhbm_gem/data/io/FileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>

namespace rg = rhbm_gem;
using MapObject = rg::MapObject;

TEST(MapFileIOTest, ReadWriteMapAndCCP4)
{
    std::array<int, 3> grid_size{ 4, 4, 4 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 1.0f, 2.0f, 3.0f };
    size_t voxel_size{ static_cast<size_t>(grid_size[0] * grid_size[1] * grid_size[2]) };
    auto values{ std::make_unique<float[]>(voxel_size) };
    for (size_t i = 0; i < voxel_size; i++)
    {
        values[i] = static_cast<float>(i);
    }
    MapObject map{ grid_size, grid_spacing, origin, std::move(values) };

    std::vector<std::string> extensions{ ".map", ".ccp4" };
    for (const auto & ext : extensions)
    {
        std::filesystem::path path{ std::filesystem::temp_directory_path() / ("test_map_io" + ext) };
        rg::WriteMap(path, map);
        auto reader{ rg::ReadMap(path) };
        ASSERT_NE(reader, nullptr);

        EXPECT_EQ(grid_size, reader->GetGridSize());
        EXPECT_EQ(grid_spacing, reader->GetGridSpacing());
        EXPECT_EQ(origin, reader->GetOrigin());

        const auto * read_values{ reader->GetMapValueArray() };
        ASSERT_NE(read_values, nullptr);
        for (size_t i = 0; i < voxel_size; i++)
        {
            EXPECT_FLOAT_EQ(static_cast<float>(i), read_values[i]);
        }
        std::filesystem::remove(path);
    }
}
