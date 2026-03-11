#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "CommandTestHelpers.hpp"

namespace {

std::vector<std::string> CollectPublicHeadersForDomain(const std::string& domain) {
    const auto include_root{command_test::ProjectRootPath() / "include" / "rhbm_gem"};
    const auto domain_root{include_root / domain};
    std::vector<std::string> headers;

    for (const auto& entry : std::filesystem::recursive_directory_iterator(domain_root)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        headers.push_back(std::filesystem::relative(entry.path(), include_root).generic_string());
    }
    std::sort(headers.begin(), headers.end());
    return headers;
}

} // namespace

TEST(DataPublicHeaderSurfaceTest, DataPublicHeadersMatchApprovedSurface) {
    const std::vector<std::string> expected{
        "data/dispatch/DataObjectDispatch.hpp",
        "data/io/DataIoServices.hpp",
        "data/io/DataObjectManager.hpp",
        "data/io/MapFileReader.hpp",
        "data/io/MapFileWriter.hpp",
        "data/io/ModelFileReader.hpp",
        "data/io/ModelFileWriter.hpp",
        "data/object/AtomClassifier.hpp",
        "data/object/AtomObject.hpp",
        "data/object/AtomicModelDataBlock.hpp",
        "data/object/BondClassifier.hpp",
        "data/object/BondObject.hpp",
        "data/object/ChemicalComponentEntry.hpp",
        "data/object/DataObjectBase.hpp",
        "data/object/GroupPotentialEntry.hpp",
        "data/object/LocalPotentialEntry.hpp",
        "data/object/MapObject.hpp",
        "data/object/ModelObject.hpp",
        "data/object/PotentialEntryQuery.hpp"};

    EXPECT_EQ(CollectPublicHeadersForDomain("data"), expected);
}
