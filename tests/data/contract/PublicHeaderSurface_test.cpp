#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "support/PublicHeaderSurfaceTestSupport.hpp"

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

    EXPECT_EQ(contract_test_support::CollectPublicHeadersForDomain("data"), expected);
}
