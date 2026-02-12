#include <gtest/gtest.h>
#include <filesystem>
#include <unordered_map>

#include "DataObjectManager.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"

TEST(DataObjectManagerTest, LoadCifFile)
{
    DataObjectManager manager;
    auto cif_path{ std::filesystem::path(__FILE__).parent_path() / "data/test_model.cif" };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    EXPECT_EQ(model->GetNumberOfAtom(), 1);
}

TEST(DataObjectManagerTest, LoadCifFileWithKeyValueEntityMetadata)
{
    DataObjectManager manager;
    auto cif_path{
        std::filesystem::path(__FILE__).parent_path() / "data/test_model_keyvalue_entity.cif"
    };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    ASSERT_EQ(model->GetNumberOfAtom(), 1);

    for (const auto & atom : model->GetAtomList()) atom->SetSelectedFlag(true);
    model->FilterAtomFromSymmetry(false);
    model->Update();
    EXPECT_EQ(model->GetNumberOfSelectedAtom(), 1);
}

TEST(DataObjectManagerTest, FilterAtomFromSymmetrySkipsWhenChainMapMissing)
{
    DataObjectManager manager;
    auto cif_path{ std::filesystem::path(__FILE__).parent_path() / "data/test_model.cif" };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    ASSERT_EQ(model->GetNumberOfAtom(), 1);

    for (const auto & atom : model->GetAtomList()) atom->SetSelectedFlag(true);
    model->SetChainIDListMap(std::unordered_map<std::string, std::vector<std::string>>{});
    ASSERT_NO_THROW(model->FilterAtomFromSymmetry(false));
    model->Update();
    EXPECT_EQ(model->GetNumberOfSelectedAtom(), 1);
}
