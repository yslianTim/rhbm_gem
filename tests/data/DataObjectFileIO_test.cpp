#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include "io/file/MapAxisOrderHelper.hpp"
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>

namespace rg = rhbm_gem;

namespace {

struct FileIoFailureCase
{
    const char * name;
    std::function<void()> run;
};

} // namespace

TEST(DataObjectFileIOTest, MapReadWriteFormatMatrix)
{
    const command_test::ScopedTempDir temp_dir{ "data_runtime_map_formats" };
    std::array<int, 3> grid_size{ 4, 4, 4 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 1.0f, 2.0f, 3.0f };
    const size_t voxel_size{
        static_cast<size_t>(grid_size[0] * grid_size[1] * grid_size[2]) };
    auto values{ std::make_unique<float[]>(voxel_size) };
    for (size_t i = 0; i < voxel_size; ++i)
    {
        values[i] = static_cast<float>(i);
    }
    rg::MapObject map{ grid_size, grid_spacing, origin, std::move(values) };

    for (const std::string & extension : { ".map", ".ccp4" })
    {
        SCOPED_TRACE(extension);
        const auto path{ temp_dir.path() / ("runtime_map" + extension) };
        ASSERT_NO_THROW(rg::WriteMap(path, map));
        auto loaded_map{ rg::ReadMap(path) };
        ASSERT_NE(loaded_map, nullptr);
        EXPECT_EQ(loaded_map->GetGridSize(), grid_size);
        EXPECT_EQ(loaded_map->GetGridSpacing(), grid_spacing);
        EXPECT_EQ(loaded_map->GetOrigin(), origin);

        const auto * loaded_values{ loaded_map->GetMapValueArray() };
        ASSERT_NE(loaded_values, nullptr);
        for (size_t i = 0; i < voxel_size; ++i)
        {
            EXPECT_FLOAT_EQ(loaded_values[i], static_cast<float>(i));
        }
    }
}

TEST(DataObjectFileIOTest, MapAxisOrderCanonicalOrderReturnsOriginalBufferAndPreservesValues)
{
    std::array<int, 3> array_size{ 2, 2, 2 };
    std::array<int, 3> axis_order{ 1, 2, 3 };
    auto raw_data{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        raw_data[i] = static_cast<float>(i + 1);
    }
    const auto * original_ptr{ raw_data.get() };

    auto reordered{
        rg::map_io::ReorderToCanonicalXYZ(std::move(raw_data), array_size, axis_order) };

    ASSERT_NE(reordered, nullptr);
    EXPECT_EQ(reordered.get(), original_ptr);
    for (size_t i = 0; i < 8; ++i)
    {
        EXPECT_FLOAT_EQ(reordered[i], static_cast<float>(i + 1));
    }
}

TEST(DataObjectFileIOTest, MapAxisOrderNonCanonicalOrderReordersToXYZ)
{
    std::array<int, 3> array_size{ 2, 3, 4 };
    std::array<int, 3> axis_order{ 3, 1, 2 };
    const size_t voxel_count{ 24 };
    auto raw_data{ std::make_unique<float[]>(voxel_count) };
    for (size_t i = 0; i < voxel_count; ++i)
    {
        raw_data[i] = static_cast<float>(i);
    }

    auto reordered{
        rg::map_io::ReorderToCanonicalXYZ(std::move(raw_data), array_size, axis_order) };

    ASSERT_NE(reordered, nullptr);
    EXPECT_FLOAT_EQ(reordered[0], 0.0f);
    EXPECT_FLOAT_EQ(reordered[1], 2.0f);
    EXPECT_FLOAT_EQ(reordered[2], 4.0f);
    EXPECT_FLOAT_EQ(reordered[3], 6.0f);
    EXPECT_FLOAT_EQ(reordered[4], 8.0f);
    EXPECT_FLOAT_EQ(reordered[5], 10.0f);
    EXPECT_FLOAT_EQ(reordered[6], 12.0f);
    EXPECT_FLOAT_EQ(reordered[12], 1.0f);
    EXPECT_FLOAT_EQ(reordered[23], 23.0f);
}

TEST(DataObjectFileIOTest, MapAxisOrderRejectsInvalidDimensionsAndAxisMappings)
{
    auto valid_data{ std::make_unique<float[]>(1) };
    valid_data[0] = 1.0f;
    EXPECT_THROW(
        (void)rg::map_io::ReorderToCanonicalXYZ(
            std::move(valid_data), std::array<int, 3>{ 0, 1, 1 }, std::array<int, 3>{ 1, 2, 3 }),
        std::runtime_error);

    auto invalid_axis_data{ std::make_unique<float[]>(1) };
    invalid_axis_data[0] = 1.0f;
    EXPECT_THROW(
        (void)rg::map_io::ReorderToCanonicalXYZ(
            std::move(invalid_axis_data),
            std::array<int, 3>{ 1, 1, 1 },
            std::array<int, 3>{ 1, 1, 3 }),
        std::runtime_error);
}

TEST(DataObjectFileIOTest, UppercaseExtensionsDispatchCorrectly)
{
    const command_test::ScopedTempDir temp_dir{ "data_runtime_uppercase_ext" };
    const auto source_model_path{ command_test::TestDataPath("test_model.cif") };
    const auto uppercase_model_path{
        data_test::CopyFixtureWithNewName(source_model_path, temp_dir.path() / "TEST_MODEL.MMCIF") };

    auto model{ rg::ReadModel(uppercase_model_path) };
    ASSERT_NE(model, nullptr);
    EXPECT_EQ(model->GetNumberOfAtom(), 1);

    const auto map_object{ data_test::MakeMapObject() };
    const auto uppercase_map_path{ temp_dir.path() / "TEST_MAP.MAP" };
    ASSERT_NO_THROW(rg::WriteMap(uppercase_map_path, map_object));
    auto loaded_map{ rg::ReadMap(uppercase_map_path) };
    ASSERT_NE(loaded_map, nullptr);
    EXPECT_EQ(loaded_map->GetGridSize(), map_object.GetGridSize());
}

TEST(DataObjectFileIOTest, PdbWriteRoundTripBasicFields)
{
    const command_test::ScopedTempDir temp_dir{ "data_runtime_pdb_roundtrip" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto output_path{ temp_dir.path() / "roundtrip.pdb" };

    auto source_model{ rg::ReadModel(model_path) };
    ASSERT_NO_THROW(rg::WriteModel(output_path, *source_model));
    auto roundtrip_model{ rg::ReadModel(output_path) };
    ASSERT_GT(roundtrip_model->GetNumberOfAtom(), 0);
    EXPECT_EQ(roundtrip_model->GetNumberOfAtom(), source_model->GetNumberOfAtom());
    EXPECT_EQ(
        roundtrip_model->GetAtomList().front()->GetSerialID(),
        source_model->GetAtomList().front()->GetSerialID());
    EXPECT_EQ(
        roundtrip_model->GetAtomList().front()->GetChainID(),
        source_model->GetAtomList().front()->GetChainID());
}

TEST(DataObjectFileIOTest, ModelWriteSupportMatrixAllowsPdbAndCifAndRejectsMmcif)
{
    const command_test::ScopedTempDir temp_dir{ "data_runtime_output_matrix" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto pdb_output_path{ temp_dir.path() / "supported_output.pdb" };
    const auto cif_output_path{ temp_dir.path() / "supported_output.cif" };
    const auto unsupported_output_path{ temp_dir.path() / "unsupported_output.mmcif" };

    auto model{ rg::ReadModel(model_path) };
    EXPECT_NO_THROW(rg::WriteModel(pdb_output_path, *model));
    EXPECT_NO_THROW(rg::WriteModel(cif_output_path, *model));
    EXPECT_THROW(rg::WriteModel(unsupported_output_path, *model), std::runtime_error);
}

TEST(DataObjectFileIOTest, ReadAndWriteFailuresSurfaceAsTypedRuntimeErrors)
{
    const command_test::ScopedTempDir temp_dir{ "data_runtime_file_io_failure" };
    const auto malformed_model_path{ temp_dir.path() / "bad_model.cif" };
    const auto missing_model_path{ temp_dir.path() / "missing_model.cif" };
    const auto malformed_map_path{ temp_dir.path() / "bad_map.map" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto bad_output_path{ temp_dir.path() / "missing_dir" / "output.cif" };

    {
        std::ofstream output{ malformed_model_path };
        output << "data_bad\nloop_\n_atom_site.id\n";
    }
    {
        std::ofstream output{ malformed_map_path, std::ios::binary };
        output << "bad";
    }

    auto model{ rg::ReadModel(model_path) };
    const std::vector<FileIoFailureCase> cases{
        { "MalformedModel", [&]() { (void)rg::ReadModel(malformed_model_path); } },
        { "MissingModel", [&]() { (void)rg::ReadModel(missing_model_path); } },
        { "MalformedMap", [&]() { (void)rg::ReadMap(malformed_map_path); } },
        { "UnopenableWriteTarget", [&]() { rg::WriteModel(bad_output_path, *model); } },
    };

    for (const auto & case_data : cases)
    {
        SCOPED_TRACE(case_data.name);
        EXPECT_THROW(case_data.run(), std::runtime_error);
    }
}

TEST(DataObjectFileIOTest, ReadModelDoesNotCallDisplay)
{
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    ASSERT_NO_THROW((void)rg::ReadModel(model_path));
    const auto stdout_output{ testing::internal::GetCapturedStdout() };
    const auto stderr_output{ testing::internal::GetCapturedStderr() };
    const auto combined_output{ stdout_output + stderr_output };
    EXPECT_EQ(combined_output.find("ModelObject Display"), std::string::npos);
    EXPECT_EQ(combined_output.find("MapObject Display"), std::string::npos);
}
