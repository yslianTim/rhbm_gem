#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>

#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

TEST(FilePathHelperTest, PathWithExtensionReturnsExtension)
{
    std::filesystem::path path{ "file.txt" };
    EXPECT_EQ(".txt", rhbm_gem::path_helper::GetExtension(path));
}

TEST(FilePathHelperTest, UppercaseExtensionIsNormalizedToLowercase)
{
    std::filesystem::path path{ "MODEL.MRC" };
    EXPECT_EQ(".mrc", rhbm_gem::path_helper::GetExtension(path));
}

TEST(FilePathHelperTest, PathWithoutExtensionReturnsEmptyString)
{
    std::filesystem::path path{ "file" };
    EXPECT_EQ("", rhbm_gem::path_helper::GetExtension(path));
}

TEST(FilePathHelperTest, PathWithMultipleExtensionsReturnsLastExtension)
{
    std::filesystem::path path{ "archive.tar.gz" };
    EXPECT_EQ(".gz", rhbm_gem::path_helper::GetExtension(path));
}

TEST(FilePathHelperTest, PathEndingWithDotReturnsDot)
{
    std::filesystem::path path{ "file." };
    EXPECT_EQ(".", rhbm_gem::path_helper::GetExtension(path));
}

TEST(FilePathHelperTest, HiddenFileReturnsFullNameAsExtension)
{
    std::filesystem::path path{ ".gitignore" };
    EXPECT_EQ(".gitignore", rhbm_gem::path_helper::GetExtension(path));
}

TEST(FilePathHelperTest, HiddenFileWithExtensionReturnsExtension)
{
    std::filesystem::path path{ ".bashrc.backup" };
    EXPECT_EQ(".backup", rhbm_gem::path_helper::GetExtension(path));
}

TEST(FilePathHelperTest, NestedDirectoriesReturnDirectoryWithTrailingSeparator)
{
    std::filesystem::path path{ std::filesystem::path("nested") / "dir" / "file.txt" };
    std::string expected{ (std::filesystem::path("nested") / "dir").string() };
    expected.push_back(std::filesystem::path::preferred_separator);
    EXPECT_EQ(expected, rhbm_gem::path_helper::GetDirectory(path));
}

TEST(FilePathHelperTest, PathWithoutParentDirectoryReturnsEmptyString)
{
    EXPECT_EQ("", rhbm_gem::path_helper::GetDirectory("file.txt"));
}

TEST(FilePathHelperTest, PathWithCurrentDirectoryReturnsDotSlash)
{
    EXPECT_EQ("./", rhbm_gem::path_helper::GetDirectory("./file.txt"));
}

TEST(FilePathHelperTest, EmptyPathReturnsEmptyString)
{
    EXPECT_EQ("", rhbm_gem::path_helper::GetDirectory(std::filesystem::path{}));
}

TEST(FilePathHelperTest, RootDirectoryPathReturnsSlash)
{
    EXPECT_EQ("/", rhbm_gem::path_helper::GetDirectory(std::filesystem::path("/")));
}

TEST(FilePathHelperTest, FileInRootDirectoryReturnsSlash)
{
    EXPECT_EQ("/", rhbm_gem::path_helper::GetDirectory(std::filesystem::path("/file.txt")));
}

TEST(FilePathHelperTest, DirectoryPathEndingWithSeparatorReturnsSamePath)
{
    std::string dir_path{ (std::filesystem::path("some") / "dir").string() };
    dir_path.push_back(std::filesystem::path::preferred_separator);
    EXPECT_EQ(dir_path, rhbm_gem::path_helper::GetDirectory(std::filesystem::path{ dir_path }));
}

TEST(FilePathHelperTest, DirectoryPathReturnsFileNameOnly)
{
    std::filesystem::path path{"/foo/bar/baz.txt"};
    EXPECT_EQ("baz.txt", rhbm_gem::path_helper::GetFileName(path));
}

TEST(FilePathHelperTest, WindowsPathReturnsFileNameOnly)
{
    std::filesystem::path path{R"(C:\temp\file.txt)"};
    EXPECT_EQ("file.txt", rhbm_gem::path_helper::GetFileName(path));
}

TEST(FilePathHelperTest, PathEndingWithSeparatorReturnsEmptyString)
{
    std::filesystem::path path{"/foo/bar/"};
    EXPECT_TRUE(rhbm_gem::path_helper::GetFileName(path).empty());
}

TEST(FilePathHelperTest, FileWithoutExtensionInDottedDirectoryReturnsFileName)
{
    const std::filesystem::path path{
        std::filesystem::path("dir.with.dots") / "file"
    };
    EXPECT_EQ("file", rhbm_gem::path_helper::GetFileName(path, false));
}

TEST(FilePathHelperTest, RemoveExtensionIgnoresDotsInParentDirectories)
{
    const std::filesystem::path path{
        std::filesystem::path("dir.with.dots") / "archive.tar.gz"
    };
    EXPECT_EQ("archive.tar", rhbm_gem::path_helper::GetFileName(path, false));
}

TEST(FilePathHelperTest, PathWithoutTrailingSlashAppendsSeparator)
{
    const std::filesystem::path input{ "some/path" };
    std::string expected{ "some/path" };
    expected.push_back(std::filesystem::path::preferred_separator);
    EXPECT_EQ(expected, rhbm_gem::path_helper::EnsureTrailingSlash(input));
}

TEST(FilePathHelperTest, PathWithTrailingSeparatorRemainsUnchanged)
{
    std::string path_with_sep{ "some/path" };
    path_with_sep.push_back(std::filesystem::path::preferred_separator);
    const std::filesystem::path input{ path_with_sep };
    EXPECT_EQ(path_with_sep, rhbm_gem::path_helper::EnsureTrailingSlash(input));
}

TEST(FilePathHelperTest, RootPathRemainsSingleSlash)
{
    const std::filesystem::path input{ "/" };
    EXPECT_EQ("/", rhbm_gem::path_helper::EnsureTrailingSlash(input));
}

TEST(FilePathHelperTest, PathEndingWithBackslashRemainsUnchanged)
{
    const std::filesystem::path input{ R"(C:\temp\)" };
    EXPECT_EQ(R"(C:\temp\)", rhbm_gem::path_helper::EnsureTrailingSlash(input));
}

TEST(FilePathHelperTest, EmptyPathRemainsEmpty)
{
    const std::filesystem::path input{};
    EXPECT_TRUE(rhbm_gem::path_helper::EnsureTrailingSlash(input).empty());
}

TEST(FilePathHelperTest, EnsureFileExistsReturnsTrueForExistingFile)
{
    const auto temp_file{
        std::filesystem::temp_directory_path() / "FilePathHelperExisting.tmp"
    };

    {
        std::ofstream ofs(temp_file);
        ofs << "test";
    }
    EXPECT_TRUE(rhbm_gem::path_helper::EnsureFileExists(temp_file, "Temp file"));
    std::filesystem::remove(temp_file);
}

TEST(FilePathHelperTest, EnsureFileExistsLogsErrorForMissingFile)
{
    const auto missing_file{
        std::filesystem::temp_directory_path() / "FilePathHelperMissing.tmp"
    };

    std::filesystem::remove(missing_file);
    testing::internal::CaptureStderr();
    const bool exists{ rhbm_gem::path_helper::EnsureFileExists(missing_file, "Missing file") };
    const std::string output{ testing::internal::GetCapturedStderr() };

    EXPECT_FALSE(exists);
    EXPECT_NE(output.find("[Error] Missing file does not exist: " +
              missing_file.string()), std::string::npos);
}

TEST(FilePathHelperTest, EnsureFileExistsReturnsTrueForExistingDirectory)
{
    const auto temp_dir{ std::filesystem::temp_directory_path() };

    testing::internal::CaptureStderr();
    const bool exists{ rhbm_gem::path_helper::EnsureFileExists(temp_dir, "Temp dir") };
    const std::string output{ testing::internal::GetCapturedStderr() };

    EXPECT_TRUE(exists);
    EXPECT_TRUE(output.empty());
}

TEST(FilePathHelperTest, EnsureFileExistsLogsErrorForEmptyPath)
{
    const std::filesystem::path empty_path{};

    testing::internal::CaptureStderr();
    const bool exists{ rhbm_gem::path_helper::EnsureFileExists(empty_path, "Empty path") };
    const std::string output{ testing::internal::GetCapturedStderr() };

    EXPECT_FALSE(exists);
    EXPECT_NE(output.find("[Error] Empty path does not exist: "), std::string::npos);
}
