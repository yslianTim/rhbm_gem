#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "CommandTestHelpers.hpp"

TEST(DataPublicHeaderSurfaceTest, LegacyModelObjectDaoHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "LegacyModelObjectDAO.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(DataPublicHeaderSurfaceTest, FileReaderBaseHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "FileReaderBase.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(DataPublicHeaderSurfaceTest, FileWriterBaseHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "FileWriterBase.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(DataPublicHeaderSurfaceTest, FileProcessFactoryRegistryHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "FileProcessFactoryRegistry.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(DataPublicHeaderSurfaceTest, DataObjectManagerHeaderDoesNotExposeInternalDebugAccessors)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto header_path{
        project_root / "include" / "rhbm_gem" / "data" / "io" / "DataObjectManager.hpp"
    };

    std::ifstream header_stream(header_path);
    ASSERT_TRUE(header_stream.is_open());
    const std::string header_content{
        std::istreambuf_iterator<char>(header_stream),
        std::istreambuf_iterator<char>()
    };

    EXPECT_EQ(header_content.find("GetDatabaseManager"), std::string::npos);
    EXPECT_EQ(header_content.find("GetDataObjectMap"), std::string::npos);
}

TEST(DataPublicHeaderSurfaceTest, DatabaseManagerHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "DatabaseManager.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(DataPublicHeaderSurfaceTest, ModelObjectDAOSqliteHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "ModelObjectDAOSqlite.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}
