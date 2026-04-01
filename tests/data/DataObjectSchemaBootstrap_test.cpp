#include <gtest/gtest.h>

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>

#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include "io/sqlite/SQLitePersistence.hpp"
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

namespace {

class FailingDataObject final : public rg::DataObjectBase
{
public:
    std::unique_ptr<rg::DataObjectBase> Clone() const override
    {
        return std::make_unique<FailingDataObject>(*this);
    }

    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag() const override { return m_key_tag; }

private:
    std::string m_key_tag;
};

} // namespace

TEST(DataObjectSchemaBootstrapTest, EmptyDatabaseBootstrapsNormalizedSchema)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_bootstrap" };
    const auto database_path{ temp_dir.path() / "bootstrap.sqlite" };

    rg::SQLitePersistence database_manager{ database_path };

    EXPECT_EQ(data_test::GetUserVersion(database_path), 2);
    EXPECT_TRUE(data_test::HasTable(database_path, "object_catalog"));
    EXPECT_FALSE(data_test::HasTable(database_path, "object_metadata"));
    EXPECT_TRUE(data_test::HasTable(database_path, "model_object"));
    EXPECT_TRUE(data_test::HasTable(database_path, "map_list"));
    EXPECT_EQ(data_test::GetUserVersion(database_path), 2);
}

TEST(DataObjectSchemaBootstrapTest, SaveUnsupportedTypeRollsBackCatalogMutation)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_atomic_save" };
    const auto database_path{ temp_dir.path() / "atomic.sqlite" };

    rg::SQLitePersistence database_manager{ database_path };
    FailingDataObject data_object;
    data_object.SetKeyTag("failing");

    EXPECT_THROW(database_manager.Save(data_object, "failing"), std::runtime_error);
    EXPECT_EQ(data_test::CountRows(database_path, "object_catalog"), 0);
}

TEST(DataObjectSchemaBootstrapTest, UnknownSchemaVersionThrows)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_unknown_version" };
    const auto database_path{ temp_dir.path() / "unknown.sqlite" };

    data_test::SetUserVersion(database_path, 99);
    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}

TEST(DataObjectSchemaBootstrapTest, VersionOneSchemaFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_version_one" };
    const auto database_path{ temp_dir.path() / "version_one.sqlite" };

    data_test::SetUserVersion(database_path, 1);
    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}
