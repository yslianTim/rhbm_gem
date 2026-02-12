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

TEST(DataObjectManagerTest, LoadCifFileMissingModelNumberDefaultsToOne)
{
    DataObjectManager manager;
    auto cif_path{
        std::filesystem::path(__FILE__).parent_path() / "data/test_model_missing_model_num.cif"
    };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetSerialID(), 1);
}

TEST(DataObjectManagerTest, LoadCifFileAuthOnlyAtomSiteColumns)
{
    DataObjectManager manager;
    auto cif_path{
        std::filesystem::path(__FILE__).parent_path() / "data/test_model_auth_only.cif"
    };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetComponentID(), "ALA");
    EXPECT_EQ(model->GetAtomList().front()->GetChainID(), "A");
    EXPECT_EQ(model->GetAtomList().front()->GetSequenceID(), 1);
}

TEST(DataObjectManagerTest, LoadCifFileMissingNumericUsesDefaults)
{
    DataObjectManager manager;
    auto cif_path{
        std::filesystem::path(__FILE__).parent_path() / "data/test_model_missing_numeric.cif"
    };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_FLOAT_EQ(model->GetAtomList().front()->GetOccupancy(), 1.0f);
    EXPECT_FLOAT_EQ(model->GetAtomList().front()->GetTemperature(), 0.0f);
}

TEST(DataObjectManagerTest, LoadCifFileModelTwoFallback)
{
    DataObjectManager manager;
    auto cif_path{
        std::filesystem::path(__FILE__).parent_path() / "data/test_model_model2_only.cif"
    };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetSerialID(), 1);
}

TEST(DataObjectManagerTest, LoadCifFileWithDoubleQuotedAtomId)
{
    DataObjectManager manager;
    auto cif_path{
        std::filesystem::path(__FILE__).parent_path() / "data/test_model_atom_id_double_quote.cif"
    };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetAtomID(), "CA A");
}

TEST(DataObjectManagerTest, LoadMmcifExtension)
{
    DataObjectManager manager;
    auto cif_path{ std::filesystem::path(__FILE__).parent_path() / "data/test_model.mmcif" };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    EXPECT_EQ(model->GetNumberOfAtom(), 1);
}

TEST(DataObjectManagerTest, LoadCifFileDatabaseOrderKeepsEmdb)
{
    DataObjectManager manager;
    auto cif_path{
        std::filesystem::path(__FILE__).parent_path() / "data/test_model_database_order.cif"
    };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    EXPECT_EQ(model->GetEmdID(), "EMD-1234");
}

TEST(DataObjectManagerTest, LoadCifFileAltBOnlyDoesNotThrow)
{
    DataObjectManager manager;
    auto cif_path{
        std::filesystem::path(__FILE__).parent_path() / "data/test_model_alt_b_only.cif"
    };
    ASSERT_NO_THROW(manager.ProcessFile(cif_path, "model"));
    auto model{ manager.GetTypedDataObject<ModelObject>("model") };
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetIndicator(), "B");
}
