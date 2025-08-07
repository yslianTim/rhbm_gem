#include <gtest/gtest.h>
#include <filesystem>

#include "DataObjectManager.hpp"
#include "ModelObject.hpp"

TEST(DataObjectManagerTest, LoadCifFile)
{
    DataObjectManager manager;
    auto cif_path{ std::filesystem::path(__FILE__).parent_path() / "data/test_model.cif" };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    EXPECT_EQ(model->GetNumberOfAtom(), 1);
}
