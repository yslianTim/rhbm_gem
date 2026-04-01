#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "support/PublicHeaderSurfaceTestSupport.hpp"

TEST(DataPublicSurfaceTest, DataPublicHeadersMatchApprovedSurface)
{
    const std::vector<std::string> expected{
        "data/io/DataRepository.hpp",
        "data/io/ModelMapFileIO.hpp",
        "data/object/AtomClassifier.hpp",
        "data/object/AtomObject.hpp",
        "data/object/BondClassifier.hpp",
        "data/object/BondObject.hpp",
        "data/object/ChemicalComponentEntry.hpp",
        "data/object/GaussianStatistics.hpp",
        "data/object/GroupPotentialEntry.hpp",
        "data/object/LocalPotentialEntry.hpp",
        "data/object/MapObject.hpp",
        "data/object/ModelObject.hpp" };

    EXPECT_EQ(contract_test_support::CollectPublicHeadersForDomain("data"), expected);
}
