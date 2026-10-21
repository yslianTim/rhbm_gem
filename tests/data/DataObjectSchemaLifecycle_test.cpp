#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include "io/sqlite/SQLiteWrapper.hpp"
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

namespace {

void CreateMarkedDatabase(const std::filesystem::path & path, int version)
{
    rg::SQLiteWrapper database{path};
    database.Execute("CREATE TABLE schema_marker (value INTEGER PRIMARY KEY);");
    database.Execute("INSERT INTO schema_marker(value) VALUES (1);");
    database.Execute("PRAGMA user_version = " + std::to_string(version) + ";");
}

std::string ReadFile(const std::filesystem::path & path)
{
    std::ifstream input(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

void ExpectUnsupportedVersionRejectedWithoutMutation(int version)
{
    const command_test::ScopedTempDir temp_dir{"schema_reject_" + std::to_string(version)};
    const auto path{temp_dir.path() / "schema.sqlite"};
    CreateMarkedDatabase(path, version);
    const auto before{ReadFile(path)};

    EXPECT_THROW((void)rg::DataRepository(path), std::runtime_error);
    EXPECT_EQ(ReadFile(path), before);
    EXPECT_EQ(data_test::GetUserVersion(path), version);
    EXPECT_EQ(data_test::CountRows(path, "schema_marker"), 1);
}

} // namespace

TEST(DataObjectSchemaLifecycleTest, EmptyDatabaseBootstrapsSchemaV20)
{
    const command_test::ScopedTempDir temp_dir{"schema_v20_bootstrap"};
    const auto path{temp_dir.path() / "schema.sqlite"};

    { rg::DataRepository repository{path}; }

    EXPECT_EQ(data_test::GetUserVersion(path), 20);
    for (const auto table : std::array<std::string_view, 9>{
             "model_atom", "model_bond", "model_chain_map", "model_component",
             "model_component_atom", "model_component_bond", "model_joint_result",
             "model_object", "model_stage_result"})
        EXPECT_TRUE(data_test::HasTable(path, std::string(table)));
    EXPECT_NO_THROW((void)rg::DataRepository(path));
}

TEST(DataObjectSchemaLifecycleTest, ExistingSchemaV20OpensWithoutMutation)
{
    const command_test::ScopedTempDir temp_dir{"schema_v20_open"};
    const auto path{temp_dir.path() / "schema.sqlite"};
    { rg::DataRepository repository{path}; }
    const auto before{ReadFile(path)};

    EXPECT_NO_THROW((void)rg::DataRepository(path));
    EXPECT_EQ(ReadFile(path), before);
    EXPECT_EQ(data_test::GetUserVersion(path), 20);
}

TEST(DataObjectSchemaLifecycleTest, UnsupportedVersionsAreRejectedWithoutMutation)
{
    for (const int version : {17, 18, 19, 21, 99})
        ExpectUnsupportedVersionRejectedWithoutMutation(version);
}

TEST(DataObjectSchemaLifecycleTest, NonemptyUnversionedDatabaseIsRejectedWithoutMutation)
{
    const command_test::ScopedTempDir temp_dir{"schema_unversioned"};
    const auto path{temp_dir.path() / "schema.sqlite"};
    CreateMarkedDatabase(path, 0);
    const auto before{ReadFile(path)};

    EXPECT_THROW((void)rg::DataRepository(path), std::runtime_error);
    EXPECT_EQ(ReadFile(path), before);
    EXPECT_EQ(data_test::GetUserVersion(path), 0);
    EXPECT_EQ(data_test::CountRows(path, "schema_marker"), 1);
}
