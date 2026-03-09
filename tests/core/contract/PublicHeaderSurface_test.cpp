#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "CommandTestHelpers.hpp"

TEST(PublicHeaderSurfaceTest, ExperimentalBondWorkflowHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "core" / "PotentialAnalysisBondWorkflow.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(PublicHeaderSurfaceTest, BuiltInCommandCatalogHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "core" / "BuiltInCommandCatalog.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(PublicHeaderSurfaceTest, CommandIdHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "core" / "CommandId.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(PublicHeaderSurfaceTest, CommandMetadataHeaderNoLongerExposesLegacySurfaceWrapper)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto header_path{
        project_root / "include" / "rhbm_gem" / "core" / "command" / "CommandMetadata.hpp"
    };

    std::ifstream header_stream(header_path);
    ASSERT_TRUE(header_stream.is_open());
    const std::string header_content{
        std::istreambuf_iterator<char>(header_stream),
        std::istreambuf_iterator<char>()
    };

    EXPECT_EQ(header_content.find("struct CommandSurface"), std::string::npos);
    EXPECT_EQ(header_content.find("MakeCommandSurface"), std::string::npos);
    EXPECT_NE(header_content.find("enum class CommandId"), std::string::npos);
    EXPECT_NE(header_content.find("enum class CommonOption"), std::string::npos);
}

TEST(PublicHeaderSurfaceTest, CommandBaseHeaderDoesNotExposeLegacyCallerFacingHooks)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto header_path{
        project_root / "include" / "rhbm_gem" / "core" / "command" / "CommandBase.hpp"
    };

    std::ifstream header_stream(header_path);
    ASSERT_TRUE(header_stream.is_open());
    const std::string header_content{
        std::istreambuf_iterator<char>(header_stream),
        std::istreambuf_iterator<char>()
    };
    const auto public_pos{ header_content.find("public:") };
    const auto protected_pos{ header_content.find("protected:") };
    ASSERT_NE(public_pos, std::string::npos);
    ASSERT_NE(protected_pos, std::string::npos);
    ASSERT_LT(public_pos, protected_pos);
    const std::string public_surface{
        header_content.substr(public_pos, protected_pos - public_pos)
    };

    EXPECT_EQ(public_surface.find("RegisterCLIOptions"), std::string::npos);
    EXPECT_EQ(public_surface.find("ReportValidationIssues"), std::string::npos);
    EXPECT_EQ(public_surface.find("GetOptions"), std::string::npos);
    EXPECT_EQ(public_surface.find("GetCommandId"), std::string::npos);
    EXPECT_EQ(header_content.find("ProcessTypedFile"), std::string::npos);
    EXPECT_EQ(header_content.find("OptionalProcessTypedFile"), std::string::npos);
    EXPECT_EQ(header_content.find("LoadTypedObject"), std::string::npos);
}

TEST(PublicHeaderSurfaceTest, LegacyModelObjectDaoHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "LegacyModelObjectDAO.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(PublicHeaderSurfaceTest, FileReaderBaseHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "FileReaderBase.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(PublicHeaderSurfaceTest, FileWriterBaseHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "FileWriterBase.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(PublicHeaderSurfaceTest, FileProcessFactoryRegistryHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "FileProcessFactoryRegistry.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(PublicHeaderSurfaceTest, DataObjectManagerHeaderDoesNotExposeInternalDebugAccessors)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto header_path{
        project_root / "include" / "rhbm_gem" / "core" / "command" / "DataObjectManager.hpp"
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

TEST(PublicHeaderSurfaceTest, DatabaseManagerHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "DatabaseManager.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}

TEST(PublicHeaderSurfaceTest, ModelObjectDaoSqliteHeaderIsNotPublic)
{
    const auto project_root{
        command_test::ProjectRootPath()
    };
    const auto leaked_header{
        project_root / "include" / "rhbm_gem" / "data" / "ModelObjectDaoSqlite.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}
