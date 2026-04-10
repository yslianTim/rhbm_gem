#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "support/PublicHeaderSurfaceTestSupport.hpp"

TEST(DataPublicHeadersTest, DataPublicHeadersMatchApprovedSurface)
{
    const std::vector<std::string> expected{
        "data/io/DataRepository.hpp",
        "data/io/ModelMapFileIO.hpp",
        "data/object/AtomObject.hpp",
        "data/object/BondObject.hpp",
        "data/object/ChemicalComponentEntry.hpp",
        "data/object/GaussianStatistics.hpp",
        "data/object/MapObject.hpp",
        "data/object/ModelAnalysisEditor.hpp",
        "data/object/ModelAnalysisView.hpp",
        "data/object/ModelObject.hpp" };

    EXPECT_EQ(contract_test_support::CollectPublicHeadersForDomain("data"), expected);
}
