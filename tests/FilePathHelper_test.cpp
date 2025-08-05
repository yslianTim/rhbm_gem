#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>

#include "FilePathHelper.hpp"

TEST(FilePathHelperTest, PathWithExtensionReturnsExtension)
{
    std::filesystem::path path{ "file.txt" };
    EXPECT_EQ(".txt", FilePathHelper::GetExtension(path));
}

TEST(FilePathHelperTest, PathWithoutExtensionReturnsEmptyString)
{
    std::filesystem::path path{ "file" };
    EXPECT_EQ("", FilePathHelper::GetExtension(path));
}

TEST(FilePathHelperTest, PathWithMultipleExtensionsReturnsLastExtension)
{
    std::filesystem::path path{ "archive.tar.gz" };
    EXPECT_EQ(".gz", FilePathHelper::GetExtension(path));
}

TEST(FilePathHelperTest, NestedDirectoriesReturnDirectoryWithTrailingSeparator)
{
    std::filesystem::path path{ std::filesystem::path("nested") / "dir" / "file.txt" };
    std::string expected{ (std::filesystem::path("nested") / "dir").string() };
    expected.push_back(std::filesystem::path::preferred_separator);
    EXPECT_EQ(expected, FilePathHelper::GetDirectory(path));
}

TEST(FilePathHelperTest, PathWithoutParentDirectoryReturnsEmptyString)
{
    EXPECT_EQ("", FilePathHelper::GetDirectory("file.txt"));
}

TEST(FilePathHelperTest, RootDirectoryPathReturnsSlash)
{
    EXPECT_EQ("/", FilePathHelper::GetDirectory(std::filesystem::path("/")));
}

TEST(FilePathHelperTest, DirectoryPathReturnsFileNameOnly)
{
    std::filesystem::path path{"/foo/bar/baz.txt"};
    EXPECT_EQ("baz.txt", FilePathHelper::GetFileName(path));
}

TEST(FilePathHelperTest, PathEndingWithSeparatorReturnsEmptyString)
{
    std::filesystem::path path{"/foo/bar/"};
    EXPECT_TRUE(FilePathHelper::GetFileName(path).empty());
}

TEST(FilePathHelperTest, PathWithoutTrailingSlashAppendsSeparator)
{
    const std::filesystem::path input{ "some/path" };
    std::string expected{ "some/path" };
    expected.push_back(std::filesystem::path::preferred_separator);
    EXPECT_EQ(expected, FilePathHelper::EnsureTrailingSlash(input));
}

TEST(FilePathHelperTest, PathWithTrailingSeparatorRemainsUnchanged)
{
    std::string path_with_sep{ "some/path" };
    path_with_sep.push_back(std::filesystem::path::preferred_separator);
    const std::filesystem::path input{ path_with_sep };
    EXPECT_EQ(path_with_sep, FilePathHelper::EnsureTrailingSlash(input));
}

TEST(FilePathHelperTest, EmptyPathRemainsEmpty)
{
    const std::filesystem::path input{};
    EXPECT_TRUE(FilePathHelper::EnsureTrailingSlash(input).empty());
}

TEST(FilePathHelperTest, EnsureFileExistsReturnsTrueForExistingFile)
{
    const auto temp_file{
        std::filesystem::temp_directory_path() /
        "FilePathHelperExisting.tmp"
    };

    {
        std::ofstream ofs(temp_file);
        ofs << "test";
    }

    EXPECT_TRUE(FilePathHelper::EnsureFileExists(temp_file, "Temp file"));

    std::filesystem::remove(temp_file);
}

TEST(FilePathHelperTest, EnsureFileExistsLogsErrorForMissingFile)
{
    const auto missing_file{
        std::filesystem::temp_directory_path() /
        "FilePathHelperMissing.tmp"
    };

    std::filesystem::remove(missing_file);

    testing::internal::CaptureStderr();
    const bool exists{
        FilePathHelper::EnsureFileExists(missing_file, "Missing file")
    };
    const std::string output{ testing::internal::GetCapturedStderr() };

    EXPECT_FALSE(exists);
    EXPECT_NE(output.find("[Error] Missing file does not exist: " +
                      missing_file.string()),
              std::string::npos);
}
