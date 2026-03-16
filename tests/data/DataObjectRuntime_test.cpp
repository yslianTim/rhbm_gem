#include <gtest/gtest.h>

#include <array>
#include <chrono>
#include <condition_variable>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <rhbm_gem/core/command/MapSampling.hpp>
#include <rhbm_gem/core/painter/AtomPainter.hpp>
#include <rhbm_gem/core/painter/ComparisonPainter.hpp>
#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/io/FileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/DataObjectDispatch.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include <rhbm_gem/utils/math/SamplerBase.hpp>
#include "CommandTestHelpers.hpp"
#include "command/CommandDataSupport.hpp"
#include "support/PublicHeaderSurfaceTestSupport.hpp"

namespace rg = rhbm_gem;

namespace {

class BlockingModelCallback {
  public:
    void operator()(rg::DataObjectBase& data_object) {
        auto* model{rg::AsModelObject(data_object)};
        if (model == nullptr) {
            return;
        }

        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_started = true;
            m_key = model->GetKeyTag();
        }
        m_cv_started.notify_all();

        std::unique_lock<std::mutex> lock(m_mutex);
        m_cv_release.wait(lock, [this] { return m_released; });
    }

    bool WaitStarted(std::chrono::milliseconds timeout) {
        std::unique_lock<std::mutex> lock(m_mutex);
        return m_cv_started.wait_for(lock, timeout, [this] { return m_started; });
    }

    void Release() {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_released = true;
        }
        m_cv_release.notify_all();
    }

    std::string Key() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_key;
    }

  private:
    mutable std::mutex m_mutex;
    std::condition_variable m_cv_started;
    std::condition_variable m_cv_release;
    bool m_started{false};
    bool m_released{false};
    std::string m_key;
};

class SinglePointSampler : public ::SamplerBase {
  public:
    void Print() const override {}

    std::vector<std::tuple<float, std::array<float, 3>>> GenerateSamplingPoints(
        const std::array<float, 3>& position,
        const std::array<float, 3>& axis_vector) const override {
        (void)axis_vector;
        return {std::make_tuple(0.0f, position)};
    }

    unsigned int GetSamplingSize() const override { return m_sampling_size; }
    void SetSamplingSize(unsigned int value) override { m_sampling_size = value; }

  private:
    unsigned int m_sampling_size{1};
};

std::shared_ptr<rg::ModelObject> LoadModelFixture(
    rg::DataObjectManager& manager,
    const std::filesystem::path& model_path,
    const std::string& key_tag = "model") {
    manager.ProcessFile(model_path, key_tag);
    return manager.GetTypedDataObject<rg::ModelObject>(key_tag);
}

std::unique_ptr<rg::ModelObject> MakeModelWithBond() {
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.reserve(2);

    auto atom_1{std::make_unique<rg::AtomObject>()};
    atom_1->SetSerialID(1);
    atom_1->SetPosition(0.0f, 0.0f, 0.0f);

    auto atom_2{std::make_unique<rg::AtomObject>()};
    atom_2->SetSerialID(2);
    atom_2->SetPosition(1.0f, 0.0f, 0.0f);

    atom_list.emplace_back(std::move(atom_1));
    atom_list.emplace_back(std::move(atom_2));
    auto model{std::make_unique<rg::ModelObject>(std::move(atom_list))};

    std::vector<std::unique_ptr<rg::BondObject>> bond_list;
    bond_list.emplace_back(std::make_unique<rg::BondObject>(
        model->GetAtomList().at(0).get(),
        model->GetAtomList().at(1).get()));
    model->SetBondList(std::move(bond_list));
    return model;
}

rg::MapObject MakeMapObject() {
    std::array<int, 3> grid_size{2, 2, 2};
    std::array<float, 3> grid_spacing{1.0f, 1.0f, 1.0f};
    std::array<float, 3> origin{0.0f, 0.0f, 0.0f};
    auto values{std::make_unique<float[]>(8)};
    for (size_t i = 0; i < 8; ++i) {
        values[i] = static_cast<float>(i + 1);
    }
    return rg::MapObject{grid_size, grid_spacing, origin, std::move(values)};
}

std::filesystem::path CopyFixtureWithNewName(
    const std::filesystem::path& source_path,
    const std::filesystem::path& output_path) {
    std::filesystem::copy_file(
        source_path,
        output_path,
        std::filesystem::copy_options::overwrite_existing);
    return output_path;
}

struct ImportRegressionCase {
    const char* name;
    std::filesystem::path path;
    std::function<void(rg::ModelObject&)> verify;
};

} // namespace

TEST(DataObjectPublicSurfaceTest, DataPublicHeadersMatchApprovedSurface) {
    const std::vector<std::string> expected{
        "data/io/DataObjectManager.hpp",
        "data/io/FileIO.hpp",
        "data/object/AtomClassifier.hpp",
        "data/object/AtomObject.hpp",
        "data/object/BondClassifier.hpp",
        "data/object/BondObject.hpp",
        "data/object/ChemicalComponentEntry.hpp",
        "data/object/DataObjectBase.hpp",
        "data/object/DataObjectDispatch.hpp",
        "data/object/GroupPotentialEntry.hpp",
        "data/object/LocalPotentialEntry.hpp",
        "data/object/MapObject.hpp",
        "data/object/ModelObject.hpp",
        "data/object/PotentialEntryQuery.hpp"};

    EXPECT_EQ(contract_test_support::CollectPublicHeadersForDomain("data"), expected);
}

TEST(DataObjectFileIOTest, MapReadWriteFormatMatrix) {
    const command_test::ScopedTempDir temp_dir{"data_runtime_map_formats"};
    std::array<int, 3> grid_size{4, 4, 4};
    std::array<float, 3> grid_spacing{1.0f, 1.0f, 1.0f};
    std::array<float, 3> origin{1.0f, 2.0f, 3.0f};
    const size_t voxel_size{
        static_cast<size_t>(grid_size[0] * grid_size[1] * grid_size[2])};
    auto values{std::make_unique<float[]>(voxel_size)};
    for (size_t i = 0; i < voxel_size; ++i) {
        values[i] = static_cast<float>(i);
    }
    rg::MapObject map{grid_size, grid_spacing, origin, std::move(values)};

    for (const std::string& extension : {".map", ".ccp4"}) {
        SCOPED_TRACE(extension);
        const auto path{temp_dir.path() / ("runtime_map" + extension)};
        ASSERT_NO_THROW(rg::WriteMap(path, map));
        auto loaded_map{rg::ReadMap(path)};
        ASSERT_NE(loaded_map, nullptr);
        EXPECT_EQ(loaded_map->GetGridSize(), grid_size);
        EXPECT_EQ(loaded_map->GetGridSpacing(), grid_spacing);
        EXPECT_EQ(loaded_map->GetOrigin(), origin);

        const auto* loaded_values{loaded_map->GetMapValueArray()};
        ASSERT_NE(loaded_values, nullptr);
        for (size_t i = 0; i < voxel_size; ++i) {
            EXPECT_FLOAT_EQ(loaded_values[i], static_cast<float>(i));
        }
    }
}

TEST(DataObjectFileIOTest, UppercaseExtensionsDispatchCorrectly) {
    const command_test::ScopedTempDir temp_dir{"data_runtime_uppercase_ext"};
    const auto source_model_path{command_test::TestDataPath("test_model.cif")};
    const auto uppercase_model_path{
        CopyFixtureWithNewName(source_model_path, temp_dir.path() / "TEST_MODEL.MMCIF")};

    rg::DataObjectManager manager{};
    ASSERT_NO_THROW(manager.ProcessFile(uppercase_model_path, "model"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("model")->GetNumberOfAtom(), 1);

    const auto map_object{MakeMapObject()};
    const auto uppercase_map_path{temp_dir.path() / "TEST_MAP.MAP"};
    ASSERT_NO_THROW(rg::WriteMap(uppercase_map_path, map_object));
    auto loaded_map{rg::ReadMap(uppercase_map_path)};
    ASSERT_NE(loaded_map, nullptr);
    EXPECT_EQ(loaded_map->GetGridSize(), map_object.GetGridSize());
}

TEST(DataObjectFileIOTest, PdbWriteRoundTripBasicFields) {
    const command_test::ScopedTempDir temp_dir{"data_runtime_pdb_roundtrip"};
    const auto model_path{command_test::TestDataPath("test_model.cif")};
    const auto output_path{temp_dir.path() / "roundtrip.pdb"};

    rg::DataObjectManager manager{};
    manager.ProcessFile(model_path, "source");
    auto source_model{manager.GetTypedDataObject<rg::ModelObject>("source")};
    ASSERT_NO_THROW(manager.ProduceFile(output_path, "source"));
    ASSERT_NO_THROW(manager.ProcessFile(output_path, "roundtrip"));

    auto roundtrip_model{manager.GetTypedDataObject<rg::ModelObject>("roundtrip")};
    ASSERT_GT(roundtrip_model->GetNumberOfAtom(), 0);
    EXPECT_EQ(roundtrip_model->GetNumberOfAtom(), source_model->GetNumberOfAtom());
    EXPECT_EQ(
        roundtrip_model->GetAtomList().front()->GetSerialID(),
        source_model->GetAtomList().front()->GetSerialID());
    EXPECT_EQ(
        roundtrip_model->GetAtomList().front()->GetChainID(),
        source_model->GetAtomList().front()->GetChainID());
}

TEST(DataObjectFileIOTest, ModelWriteSupportMatrixAllowsPdbAndCifAndRejectsMmcif) {
    const command_test::ScopedTempDir temp_dir{"data_runtime_output_matrix"};
    const auto model_path{command_test::TestDataPath("test_model.cif")};
    const auto pdb_output_path{temp_dir.path() / "supported_output.pdb"};
    const auto cif_output_path{temp_dir.path() / "supported_output.cif"};
    const auto unsupported_output_path{temp_dir.path() / "unsupported_output.mmcif"};

    rg::DataObjectManager manager{};
    manager.ProcessFile(model_path, "model");
    EXPECT_NO_THROW(manager.ProduceFile(pdb_output_path, "model"));
    EXPECT_NO_THROW(manager.ProduceFile(cif_output_path, "model"));
    EXPECT_THROW(manager.ProduceFile(unsupported_output_path, "model"), std::runtime_error);
}

TEST(DataObjectFileIOTest, ProcessFileThrowsOnMalformedModelInput) {
    const command_test::ScopedTempDir temp_dir{"data_runtime_bad_model"};
    const auto malformed_path{temp_dir.path() / "bad_model.cif"};
    {
        std::ofstream output{malformed_path};
        output << "data_bad\nloop_\n_atom_site.id\n";
    }

    rg::DataObjectManager manager{};
    EXPECT_THROW(manager.ProcessFile(malformed_path, "broken"), std::runtime_error);
}

TEST(DataObjectFileIOTest, ProcessFileThrowsOnMalformedMapInput) {
    const command_test::ScopedTempDir temp_dir{"data_runtime_bad_map"};
    const auto malformed_path{temp_dir.path() / "bad_map.map"};
    {
        std::ofstream output{malformed_path, std::ios::binary};
        output << "bad";
    }

    rg::DataObjectManager manager{};
    EXPECT_THROW(manager.ProcessFile(malformed_path, "broken_map"), std::runtime_error);
}

TEST(DataObjectFileIOTest, ProduceFileThrowsWhenWriterCannotOpenTarget) {
    const command_test::ScopedTempDir temp_dir{"data_runtime_bad_output_target"};
    const auto model_path{command_test::TestDataPath("test_model.cif")};
    const auto output_path{temp_dir.path() / "missing_dir" / "output.cif"};

    rg::DataObjectManager manager{};
    manager.ProcessFile(model_path, "model");
    EXPECT_THROW(manager.ProduceFile(output_path, "model"), std::runtime_error);
}

TEST(DataObjectFileIOTest, FunctionFileIoThrowsWhenReadFails) {
    const command_test::ScopedTempDir temp_dir{"data_runtime_file_io_failure"};
    const auto missing_model_path{temp_dir.path() / "missing_model.cif"};
    const auto malformed_map_path{temp_dir.path() / "bad_map.map"};
    {
        std::ofstream output{malformed_map_path, std::ios::binary};
        output << "bad";
    }

    EXPECT_THROW((void)rg::ReadModel(missing_model_path), std::runtime_error);
    EXPECT_THROW((void)rg::ReadMap(malformed_map_path), std::runtime_error);
}

TEST(DataObjectFileIOTest, WriteDataObjectRejectsUnsupportedTopLevelObjectType) {
    rg::AtomObject atom;
    EXPECT_THROW(rg::WriteDataObject("unused.pdb", atom), std::runtime_error);
    EXPECT_THROW(rg::WriteDataObject("unused.map", atom), std::runtime_error);
}

TEST(DataObjectImportRegressionTest, CifEdgeCaseMatrix) {
    const command_test::ScopedTempDir temp_dir{"data_runtime_mmcif_fixture"};
    const auto mmcif_path{
        CopyFixtureWithNewName(
            command_test::TestDataPath("test_model.cif"),
            temp_dir.path() / "test_model_runtime.mmcif")};

    const std::vector<ImportRegressionCase> cases{
        {"CifExtensionLoadsOneAtom",
         command_test::TestDataPath("test_model.cif"),
         [](rg::ModelObject& model) {
             EXPECT_EQ(model.GetNumberOfAtom(), 1);
         }},
        {"MmcifExtensionLoadsOneAtom",
         mmcif_path,
         [](rg::ModelObject& model) {
             EXPECT_EQ(model.GetNumberOfAtom(), 1);
         }},
        {"MissingModelNumberDefaultsSerialId",
         command_test::TestDataPath("test_model_missing_model_num.cif"),
         [](rg::ModelObject& model) {
             ASSERT_EQ(model.GetNumberOfAtom(), 1);
             EXPECT_EQ(model.GetAtomList().front()->GetSerialID(), 1);
         }},
        {"ModelTwoUsesSerialFallback",
         command_test::TestDataPath("test_model_model2_only.cif"),
         [](rg::ModelObject& model) {
             ASSERT_EQ(model.GetNumberOfAtom(), 1);
             EXPECT_EQ(model.GetAtomList().front()->GetSerialID(), 1);
         }},
        {"DoubleQuotedAtomIdPreserved",
         command_test::TestDataPath("test_model_atom_id_double_quote.cif"),
         [](rg::ModelObject& model) {
             ASSERT_EQ(model.GetNumberOfAtom(), 1);
             EXPECT_EQ(model.GetAtomList().front()->GetAtomID(), "CA A");
         }},
        {"LoopMultilineQuotedTokenPreserved",
         command_test::TestDataPath("test_model_loop_multiline.cif"),
         [](rg::ModelObject& model) {
             ASSERT_EQ(model.GetNumberOfAtom(), 1);
             EXPECT_EQ(model.GetAtomList().front()->GetAtomID(), "CA A");
         }},
        {"KeyValueEntityMetadataSupportsSymmetrySelection",
         command_test::TestDataPath("test_model_keyvalue_entity.cif"),
         [](rg::ModelObject& model) {
             ASSERT_EQ(model.GetNumberOfAtom(), 1);
             for (const auto& atom : model.GetAtomList()) {
                 atom->SetSelectedFlag(true);
             }
             model.FilterAtomFromSymmetry(false);
             model.Update();
             EXPECT_EQ(model.GetNumberOfSelectedAtom(), 1);
         }},
        {"MissingChainMapSkipsSymmetryFiltering",
         command_test::TestDataPath("test_model.cif"),
         [](rg::ModelObject& model) {
             ASSERT_EQ(model.GetNumberOfAtom(), 1);
             for (const auto& atom : model.GetAtomList()) {
                 atom->SetSelectedFlag(true);
             }
             model.SetChainIDListMap(
                 std::unordered_map<std::string, std::vector<std::string>>{});
             model.FilterAtomFromSymmetry(false);
             model.Update();
             EXPECT_EQ(model.GetNumberOfSelectedAtom(), 1);
         }},
        {"AuthOnlyColumnsPopulateFallbackFields",
         command_test::TestDataPath("test_model_auth_only.cif"),
         [](rg::ModelObject& model) {
             ASSERT_EQ(model.GetNumberOfAtom(), 1);
             EXPECT_EQ(model.GetAtomList().front()->GetComponentID(), "ALA");
             EXPECT_EQ(model.GetAtomList().front()->GetChainID(), "A");
             EXPECT_EQ(model.GetAtomList().front()->GetSequenceID(), 1);
         }},
        {"MissingNumericFieldsUseDefaults",
         command_test::TestDataPath("test_model_missing_numeric.cif"),
         [](rg::ModelObject& model) {
             ASSERT_EQ(model.GetNumberOfAtom(), 1);
             EXPECT_FLOAT_EQ(model.GetAtomList().front()->GetOccupancy(), 1.0f);
             EXPECT_FLOAT_EQ(model.GetAtomList().front()->GetTemperature(), 0.0f);
         }},
        {"DatabaseOrderKeepsEmdbIdentifier",
         command_test::TestDataPath("test_model_database_order.cif"),
         [](rg::ModelObject& model) {
             EXPECT_EQ(model.GetEmdID(), "EMD-1234");
         }},
        {"AltBIndicatorPreserved",
         command_test::TestDataPath("test_model_alt_b_only.cif"),
         [](rg::ModelObject& model) {
             ASSERT_EQ(model.GetNumberOfAtom(), 1);
             EXPECT_EQ(model.GetAtomList().front()->GetIndicator(), "B");
         }},
        {"InvalidSecondaryRangeDoesNotDropAtoms",
         command_test::TestDataPath("test_model_invalid_secondary_range.cif"),
         [](rg::ModelObject& model) {
             EXPECT_EQ(model.GetNumberOfAtom(), 1);
         }},
        {"StructConnImportBuildsBondAndBondKeySystem",
         command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif"),
         [](rg::ModelObject& model) {
             ASSERT_EQ(model.GetNumberOfAtom(), 2);
             ASSERT_GE(model.GetNumberOfBond(), 1);
             ASSERT_NE(model.GetBondKeySystemPtr(), nullptr);
             EXPECT_NE(
                 model.GetBondKeySystemPtr()->GetBondId(
                     model.GetBondList().front()->GetBondKey()),
                 "UNK");
         }}};

    for (const auto& case_data : cases) {
        SCOPED_TRACE(case_data.name);
        rg::DataObjectManager manager{};
        auto model{LoadModelFixture(manager, case_data.path)};
        case_data.verify(*model);
    }
}

TEST(DataObjectManagerIterationTest, ForEachDataObjectUsesDeterministicOrderByDefault) {
    rg::DataObjectManager manager{};
    const auto model_path{command_test::TestDataPath("test_model.cif")};
    manager.ProcessFile(model_path, "b_model");
    manager.ProcessFile(model_path, "a_model");

    rg::DataObjectManager::IterateOptions options;
    std::vector<std::string> keys;
    manager.ForEachDataObject(
        [&keys](const rg::DataObjectBase& data_object) {
            if (const auto* model{rg::AsModelObject(data_object)}) {
                keys.push_back(model->GetKeyTag());
            }
        },
        {},
        options);

    EXPECT_EQ(keys, (std::vector<std::string>{"a_model", "b_model"}));
}

TEST(DataObjectManagerIterationTest, ForEachDataObjectPreservesExplicitKeyOrder) {
    rg::DataObjectManager manager{};
    const auto model_path{command_test::TestDataPath("test_model.cif")};
    manager.ProcessFile(model_path, "b_model");
    manager.ProcessFile(model_path, "a_model");

    std::vector<std::string> keys;
    manager.ForEachDataObject(
        [&keys](const rg::DataObjectBase& data_object) {
            if (const auto* model{rg::AsModelObject(data_object)}) {
                keys.push_back(model->GetKeyTag());
            }
        },
        {"b_model", "a_model"});

    EXPECT_EQ(keys, (std::vector<std::string>{"b_model", "a_model"}));
}

TEST(DataObjectManagerIterationTest, ConstOverloadWorks) {
    rg::DataObjectManager manager{};
    const auto model_path{command_test::TestDataPath("test_model.cif")};
    manager.ProcessFile(model_path, "b_model");
    manager.ProcessFile(model_path, "a_model");

    const rg::DataObjectManager& const_manager{manager};
    std::vector<std::string> keys;
    const_manager.ForEachDataObject(
        [&keys](const rg::DataObjectBase& data_object) {
            if (const auto* model{rg::AsModelObject(data_object)}) {
                keys.push_back(model->GetKeyTag());
            }
        });

    EXPECT_EQ(keys, (std::vector<std::string>{"a_model", "b_model"}));
}

TEST(DataObjectManagerIterationTest, SnapshotAllowsClearDuringCallback) {
    rg::DataObjectManager manager{};
    manager.ProcessFile(command_test::TestDataPath("test_model.cif"), "model");

    BlockingModelCallback callback;
    std::exception_ptr worker_error;
    std::thread worker(
        [&] {
            try {
                manager.ForEachDataObject(std::ref(callback), {"model"});
            } catch (...) {
                worker_error = std::current_exception();
            }
        });

    ASSERT_TRUE(callback.WaitStarted(std::chrono::seconds(2)));
    manager.ClearDataObjects();
    callback.Release();
    worker.join();

    EXPECT_EQ(worker_error, nullptr);
    EXPECT_EQ(callback.Key(), "model");
}

TEST(DataObjectManagerIterationTest, RejectsEmptyCallback) {
    rg::DataObjectManager manager{};
    std::function<void(rg::DataObjectBase&)> mutable_callback;
    EXPECT_THROW(manager.ForEachDataObject(mutable_callback), std::runtime_error);

    const rg::DataObjectManager& const_manager{manager};
    std::function<void(const rg::DataObjectBase&)> const_callback;
    EXPECT_THROW(const_manager.ForEachDataObject(const_callback), std::runtime_error);
}

TEST(DataObjectDispatchTest, AsHelpersResolveModelAndMap) {
    auto model{MakeModelWithBond()};
    auto map{MakeMapObject()};
    rg::AtomObject atom;

    EXPECT_EQ(rg::AsModelObject(*model), model.get());
    EXPECT_EQ(rg::AsMapObject(*model), nullptr);
    EXPECT_EQ(rg::AsMapObject(map), &map);
    EXPECT_EQ(rg::AsModelObject(map), nullptr);
    EXPECT_EQ(rg::AsModelObject(atom), nullptr);
    EXPECT_EQ(rg::AsMapObject(atom), nullptr);

    const rg::DataObjectBase& const_model_ref{*model};
    const rg::DataObjectBase& const_map_ref{map};
    EXPECT_EQ(rg::AsModelObject(const_model_ref), model.get());
    EXPECT_EQ(rg::AsMapObject(const_map_ref), &map);
}

TEST(DataObjectDispatchTest, ExpectHelpersResolveModelAndMap) {
    auto model{MakeModelWithBond()};
    auto map{MakeMapObject()};

    const auto& model_ref{rg::ExpectModelObject(*model, "dispatch-test-model")};
    const auto& map_ref{rg::ExpectMapObject(map, "dispatch-test-map")};
    EXPECT_EQ(&model_ref, model.get());
    EXPECT_EQ(&map_ref, &map);
}

TEST(DataObjectDispatchTest, ExpectHelpersThrowOnUnsupportedTargetType) {
    rg::AtomObject atom;
    EXPECT_THROW((void)rg::ExpectModelObject(atom, "dispatch-test-model"), std::runtime_error);
    EXPECT_THROW((void)rg::ExpectMapObject(atom, "dispatch-test-map"), std::runtime_error);
}

TEST(DataObjectDispatchTest, CatalogTypeNameUsesStableTopLevelNames) {
    auto model{MakeModelWithBond()};
    auto map{MakeMapObject()};
    rg::AtomObject atom;

    EXPECT_EQ(rg::GetCatalogTypeName(*model), "model");
    EXPECT_EQ(rg::GetCatalogTypeName(map), "map");
    EXPECT_THROW((void)rg::GetCatalogTypeName(atom), std::runtime_error);
}

TEST(DataObjectOperationTest, SampleMapValuesReturnsExpectedPointValueAndIsDeterministic) {
    auto map{MakeMapObject()};
    SinglePointSampler sampler;
    const auto first{
        rg::SampleMapValues(map, sampler, {0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f})};
    const auto second{
        rg::SampleMapValues(map, sampler, {1.0f, 1.0f, 1.0f}, {0.0f, 1.0f, 0.0f})};
    const auto first_again{
        rg::SampleMapValues(map, sampler, {0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f})};

    ASSERT_EQ(first.size(), 1);
    ASSERT_EQ(second.size(), 1);
    ASSERT_EQ(first_again.size(), 1);
    EXPECT_FLOAT_EQ(std::get<0>(first.front()), 0.0f);
    EXPECT_FLOAT_EQ(std::get<0>(second.front()), 0.0f);
    EXPECT_FLOAT_EQ(std::get<1>(first.front()), map.GetMapValue(0, 0, 0));
    EXPECT_FLOAT_EQ(std::get<1>(second.front()), map.GetMapValue(1, 1, 1));
    EXPECT_FLOAT_EQ(std::get<1>(first_again.front()), std::get<1>(first.front()));
}

TEST(DataObjectOperationTest, AtomPainterDispatchesByTypedIngestionAndRejectsUnsupportedType) {
    rg::AtomPainter painter;
    rg::AtomObject atom;
    atom.SetSelectedFlag(true);
    atom.AddLocalPotentialEntry(std::make_unique<rg::LocalPotentialEntry>());
    painter.AddDataObject(&atom);
    painter.AddReferenceDataObject(&atom, "ref");

    auto model{MakeModelWithBond()};
    EXPECT_THROW(painter.AddDataObject(model.get()), std::runtime_error);
    EXPECT_THROW(painter.AddReferenceDataObject(model.get(), "ref"), std::runtime_error);
}

TEST(DataObjectOperationTest, ModelBasedPaintersDispatchByTypedIngestionAndRejectUnsupportedType) {
    auto model{MakeModelWithBond()};
    rg::AtomObject atom;

    rg::ModelPainter model_painter;
    model_painter.AddDataObject(model.get());
    model_painter.AddReferenceDataObject(model.get(), "ref");
    EXPECT_THROW(model_painter.AddDataObject(&atom), std::runtime_error);

    rg::GausPainter gaus_painter;
    gaus_painter.AddDataObject(model.get());
    gaus_painter.AddReferenceDataObject(model.get(), "ref");
    EXPECT_THROW(gaus_painter.AddDataObject(&atom), std::runtime_error);

    rg::ComparisonPainter comparison_painter;
    comparison_painter.AddDataObject(model.get());
    comparison_painter.AddReferenceDataObject(model.get(), "ref");
    EXPECT_THROW(comparison_painter.AddDataObject(&atom), std::runtime_error);

    rg::DemoPainter demo_painter;
    demo_painter.AddDataObject(model.get());
    demo_painter.AddReferenceDataObject(model.get(), "ref");
    EXPECT_THROW(demo_painter.AddDataObject(&atom), std::runtime_error);
}

TEST(DataObjectOperationTest, PainterNullIngestThrows) {
    rg::ModelPainter painter;
    EXPECT_THROW(painter.AddDataObject(nullptr), std::runtime_error);
    EXPECT_THROW(painter.AddReferenceDataObject(nullptr, "ref"), std::runtime_error);
}

TEST(DataObjectOperationTest, NormalizeMapObjectNormalizesMapValues) {
    auto map{MakeMapObject()};
    const auto original_value{map.GetMapValue(0)};
    const auto original_sd{map.GetMapValueSD()};
    ASSERT_GT(original_sd, 0.0f);

    rg::NormalizeMapObject(map);

    EXPECT_NEAR(map.GetMapValue(0), original_value / original_sd, 1.0e-5f);
}

TEST(DataObjectOperationTest, PrepareModelObjectSelectsAndInitializesLocalEntries) {
    auto model{MakeModelWithBond()};
    for (auto& atom : model->GetAtomList()) {
        atom->SetSelectedFlag(false);
    }
    for (auto& bond : model->GetBondList()) {
        bond->SetSelectedFlag(false);
    }
    model->Update();
    ASSERT_EQ(model->GetNumberOfSelectedAtom(), 0);
    ASSERT_EQ(model->GetNumberOfSelectedBond(), 0);

    rg::ModelPreparationOptions options;
    options.select_all_atoms = true;
    options.select_all_bonds = true;
    options.update_model = true;
    options.initialize_atom_local_entries = true;
    options.initialize_bond_local_entries = true;

    rg::PrepareModelObject(*model, options);

    EXPECT_EQ(model->GetNumberOfSelectedAtom(), model->GetNumberOfAtom());
    EXPECT_EQ(model->GetNumberOfSelectedBond(), model->GetNumberOfBond());
    for (const auto* atom : model->GetSelectedAtomList()) {
        ASSERT_NE(atom, nullptr);
        EXPECT_NE(atom->GetLocalPotentialEntry(), nullptr);
    }
    for (const auto* bond : model->GetSelectedBondList()) {
        ASSERT_NE(bond, nullptr);
        EXPECT_NE(bond->GetLocalPotentialEntry(), nullptr);
    }
}

TEST(DataObjectOperationTest, ApplyModelSelectionSelectsByAtomSelectorRules) {
    auto model{MakeModelWithBond()};
    auto& atom_list{model->GetAtomList()};
    ASSERT_EQ(atom_list.size(), 2);
    atom_list.at(0)->SetChainID("A");
    atom_list.at(1)->SetChainID("B");
    atom_list.at(0)->SetSelectedFlag(false);
    atom_list.at(1)->SetSelectedFlag(true);

    ::AtomSelector selector;
    selector.PickChainID("A");
    rg::ApplyModelSelection(*model, selector);

    EXPECT_TRUE(atom_list.at(0)->GetSelectedFlag());
    EXPECT_FALSE(atom_list.at(1)->GetSelectedFlag());
}

TEST(DataObjectOperationTest, CollectModelAtomsSupportsSelectionAndEntryFilters) {
    auto model{MakeModelWithBond()};
    auto& atoms{model->GetAtomList()};
    ASSERT_EQ(atoms.size(), 2);
    atoms[0]->SetSelectedFlag(true);
    atoms[1]->SetSelectedFlag(false);
    atoms[0]->AddLocalPotentialEntry(std::make_unique<rg::LocalPotentialEntry>());
    model->Update();

    rg::ModelAtomCollectorOptions selected_only_options;
    selected_only_options.selected_only = true;
    const auto selected_only_atoms{rg::CollectModelAtoms(*model, selected_only_options)};
    ASSERT_EQ(selected_only_atoms.size(), 1);
    EXPECT_EQ(selected_only_atoms.front(), atoms[0].get());

    rg::ModelAtomCollectorOptions require_entry_options;
    require_entry_options.selected_only = false;
    require_entry_options.require_local_potential_entry = true;
    const auto require_entry_atoms{rg::CollectModelAtoms(*model, require_entry_options)};
    ASSERT_EQ(require_entry_atoms.size(), 1);
    EXPECT_EQ(require_entry_atoms.front(), atoms[0].get());
}

TEST(DataObjectOperationTest, PrepareSimulationAtomsCollectsAtomChargeAndRange) {
    auto model{MakeModelWithBond()};
    rg::SimulationAtomPreparationOptions options;
    options.partial_charge_choice = rg::PartialCharge::NEUTRAL;
    options.include_unknown_atoms = true;

    const auto result{rg::PrepareSimulationAtoms(*model, options)};
    ASSERT_TRUE(result.has_atom);
    EXPECT_EQ(result.atom_list.size(), 2);
    ASSERT_EQ(result.atom_charge_map.size(), 2);
    EXPECT_DOUBLE_EQ(result.atom_charge_map.at(1), 0.0);
    EXPECT_DOUBLE_EQ(result.atom_charge_map.at(2), 0.0);
    EXPECT_FLOAT_EQ(result.range_minimum[0], 0.0f);
    EXPECT_FLOAT_EQ(result.range_maximum[0], 1.0f);
}

TEST(DataObjectOperationTest, BuildModelAtomBondContextBuildsSelectedContextMaps) {
    auto model{MakeModelWithBond()};
    auto& atoms{model->GetAtomList()};
    auto& bonds{model->GetBondList()};
    ASSERT_EQ(atoms.size(), 2);
    ASSERT_EQ(bonds.size(), 1);
    atoms[0]->SetSelectedFlag(true);
    atoms[1]->SetSelectedFlag(true);
    bonds[0]->SetSelectedFlag(true);
    model->Update();

    const auto context{rg::BuildModelAtomBondContext(*model)};
    ASSERT_EQ(context.atom_map.size(), 2);
    ASSERT_EQ(context.bond_map.size(), 2);
    EXPECT_EQ(context.atom_map.at(1), atoms[0].get());
    EXPECT_EQ(context.atom_map.at(2), atoms[1].get());
    EXPECT_EQ(context.bond_map.at(1).size(), 1);
    EXPECT_EQ(context.bond_map.at(2).size(), 1);
}
