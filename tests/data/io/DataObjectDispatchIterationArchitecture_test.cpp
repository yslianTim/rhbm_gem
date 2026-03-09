#include <gtest/gtest.h>

#include <array>
#include <chrono>
#include <condition_variable>
#include <exception>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/core/painter/AtomPainter.hpp>
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/core/painter/ComparisonPainter.hpp>
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/data/dispatch/DataObjectDispatch.hpp>
#include <rhbm_gem/core/command/DataObjectManager.hpp>
#include "workflow/DataObjectWorkflowOps.hpp"
#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include "support/DataInternalTestSeam.hpp"
#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/core/command/MapSampling.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/utils/math/SamplerBase.hpp>

namespace rg = rhbm_gem;

namespace {

class BlockingModelCallback
{
public:
    void operator()(rg::DataObjectBase & data_object)
    {
        auto * model{ rg::AsModelObject(data_object) };
        if (model == nullptr) return;

        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_started = true;
            m_key = model->GetKeyTag();
        }
        m_cv_started.notify_all();

        std::unique_lock<std::mutex> lock(m_mutex);
        m_cv_release.wait(lock, [this] { return m_released; });
    }

    bool WaitStarted(std::chrono::milliseconds timeout)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        return m_cv_started.wait_for(lock, timeout, [this] { return m_started; });
    }

    void Release()
    {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_released = true;
        }
        m_cv_release.notify_all();
    }

    std::string Key() const
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_key;
    }

private:
    mutable std::mutex m_mutex;
    std::condition_variable m_cv_started;
    std::condition_variable m_cv_release;
    bool m_started{ false };
    bool m_released{ false };
    std::string m_key;
};

class SinglePointSampler : public ::SamplerBase
{
public:
    void Print() const override {}

    std::vector<std::tuple<float, std::array<float, 3>>> GenerateSamplingPoints(
        const std::array<float, 3> & position,
        const std::array<float, 3> & axis_vector) const override
    {
        (void)axis_vector;
        return { std::make_tuple(0.0f, position) };
    }

    unsigned int GetSamplingSize() const override { return m_sampling_size; }
    void SetSamplingSize(unsigned int value) override { m_sampling_size = value; }

private:
    unsigned int m_sampling_size{ 1 };
};

std::unique_ptr<rg::ModelObject> MakeModelWithBond()
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.reserve(2);

    auto atom_1{ std::make_unique<rg::AtomObject>() };
    atom_1->SetSerialID(1);
    atom_1->SetPosition(0.0f, 0.0f, 0.0f);

    auto atom_2{ std::make_unique<rg::AtomObject>() };
    atom_2->SetSerialID(2);
    atom_2->SetPosition(1.0f, 0.0f, 0.0f);

    atom_list.emplace_back(std::move(atom_1));
    atom_list.emplace_back(std::move(atom_2));
    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };

    std::vector<std::unique_ptr<rg::BondObject>> bond_list;
    bond_list.emplace_back(std::make_unique<rg::BondObject>(
        model->GetAtomList().at(0).get(), model->GetAtomList().at(1).get()));
    model->SetBondList(std::move(bond_list));
    return model;
}

rg::MapObject MakeMapObject()
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; i++) values[i] = static_cast<float>(i + 1);
    return rg::MapObject{ grid_size, grid_spacing, origin, std::move(values) };
}

} // namespace

TEST(DataObjectDispatchIterationArchitectureTest, ManagerForEachDataObjectDeterministicOrderByDefault)
{
    rg::DataObjectManager manager;
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    manager.ProcessFile(model_path, "b_model");
    manager.ProcessFile(model_path, "a_model");

    rg::DataObjectManager::IterateOptions options;
    std::vector<std::string> keys;
    manager.ForEachDataObject(
        [&keys](const rg::DataObjectBase & data_object)
        {
            if (const auto * model{ rg::AsModelObject(data_object) })
            {
                keys.push_back(model->GetKeyTag());
            }
        },
        {},
        options);

    EXPECT_EQ(keys, (std::vector<std::string>{ "a_model", "b_model" }));
}

TEST(DataObjectDispatchIterationArchitectureTest, ManagerForEachDataObjectPreservesInputKeyOrder)
{
    rg::DataObjectManager manager;
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    manager.ProcessFile(model_path, "b_model");
    manager.ProcessFile(model_path, "a_model");

    std::vector<std::string> keys;
    manager.ForEachDataObject(
        [&keys](const rg::DataObjectBase & data_object)
        {
            if (const auto * model{ rg::AsModelObject(data_object) })
            {
                keys.push_back(model->GetKeyTag());
            }
        },
        { "b_model", "a_model" });

    EXPECT_EQ(keys, (std::vector<std::string>{ "b_model", "a_model" }));
}

TEST(DataObjectDispatchIterationArchitectureTest, ConstManagerForEachDataObjectWorks)
{
    rg::DataObjectManager manager;
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    manager.ProcessFile(model_path, "b_model");
    manager.ProcessFile(model_path, "a_model");

    const rg::DataObjectManager & const_manager{ manager };
    std::vector<std::string> keys;
    const_manager.ForEachDataObject(
        [&keys](const rg::DataObjectBase & data_object)
        {
            if (const auto * model{ rg::AsModelObject(data_object) })
            {
                keys.push_back(model->GetKeyTag());
            }
        });

    EXPECT_EQ(keys, (std::vector<std::string>{ "a_model", "b_model" }));
}

TEST(DataObjectDispatchIterationArchitectureTest, ManagerForEachDataObjectUsesSnapshotSoClearCanRunConcurrently)
{
    rg::DataObjectManager manager;
    manager.ProcessFile(command_test::TestDataPath("test_model.cif"), "model");

    BlockingModelCallback callback;
    std::exception_ptr worker_error;
    std::thread worker(
        [&]
        {
            try
            {
                manager.ForEachDataObject(std::ref(callback), { "model" });
            }
            catch (...)
            {
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

TEST(DataObjectDispatchIterationArchitectureTest, ManagerForEachDataObjectRejectsEmptyCallback)
{
    rg::DataObjectManager manager;
    std::function<void(rg::DataObjectBase &)> mutable_callback;
    EXPECT_THROW(manager.ForEachDataObject(mutable_callback), std::runtime_error);

    const rg::DataObjectManager & const_manager{ manager };
    std::function<void(const rg::DataObjectBase &)> const_callback;
    EXPECT_THROW(const_manager.ForEachDataObject(const_callback), std::runtime_error);
}

TEST(DataObjectDispatchIterationArchitectureTest, SampleMapValuesProducesSamplingOutput)
{
    auto map{ MakeMapObject() };
    SinglePointSampler sampler;
    const auto sampling_data{
        rg::SampleMapValues(map, sampler, { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f })
    };

    ASSERT_EQ(sampling_data.size(), 1);
    EXPECT_FLOAT_EQ(std::get<0>(sampling_data.front()), 0.0f);
}

TEST(DataObjectDispatchIterationArchitectureTest, SampleMapValuesIsStatelessAcrossCalls)
{
    auto map{ MakeMapObject() };
    SinglePointSampler sampler;
    const auto first{
        rg::SampleMapValues(map, sampler, { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f })
    };
    const auto second{
        rg::SampleMapValues(map, sampler, { 1.0f, 1.0f, 1.0f }, { 0.0f, 1.0f, 0.0f })
    };

    EXPECT_EQ(first.size(), 1);
    EXPECT_EQ(second.size(), 1);
}

TEST(DataObjectDispatchIterationArchitectureTest, AtomPainterDispatchesByTypedIngestionAndRejectsUnsupportedType)
{
    rg::AtomPainter painter;
    rg::AtomObject atom;
    atom.SetSelectedFlag(true);
    atom.AddLocalPotentialEntry(std::make_unique<rg::LocalPotentialEntry>());
    EXPECT_NO_THROW(painter.AddDataObject(&atom));
    EXPECT_NO_THROW(painter.AddReferenceDataObject(&atom, "ref"));

    auto model{ MakeModelWithBond() };
    EXPECT_THROW(painter.AddDataObject(model.get()), std::runtime_error);
    EXPECT_THROW(painter.AddReferenceDataObject(model.get(), "ref"), std::runtime_error);
}

TEST(DataObjectDispatchIterationArchitectureTest, ModelBasedPaintersDispatchByTypedIngestionAndRejectUnsupportedType)
{
    auto model{ MakeModelWithBond() };
    rg::AtomObject atom;

    rg::ModelPainter model_painter;
    EXPECT_NO_THROW(model_painter.AddDataObject(model.get()));
    EXPECT_NO_THROW(model_painter.AddReferenceDataObject(model.get(), "ref"));
    EXPECT_THROW(model_painter.AddDataObject(&atom), std::runtime_error);

    rg::GausPainter gaus_painter;
    EXPECT_NO_THROW(gaus_painter.AddDataObject(model.get()));
    EXPECT_NO_THROW(gaus_painter.AddReferenceDataObject(model.get(), "ref"));
    EXPECT_THROW(gaus_painter.AddDataObject(&atom), std::runtime_error);

    rg::ComparisonPainter comparison_painter;
    EXPECT_NO_THROW(comparison_painter.AddDataObject(model.get()));
    EXPECT_NO_THROW(comparison_painter.AddReferenceDataObject(model.get(), "ref"));
    EXPECT_THROW(comparison_painter.AddDataObject(&atom), std::runtime_error);

    rg::DemoPainter demo_painter;
    EXPECT_NO_THROW(demo_painter.AddDataObject(model.get()));
    EXPECT_NO_THROW(demo_painter.AddReferenceDataObject(model.get(), "ref"));
    EXPECT_THROW(demo_painter.AddDataObject(&atom), std::runtime_error);
}

TEST(DataObjectDispatchIterationArchitectureTest, PainterNullIngestUsesModeSpecificErrorMessage)
{
    rg::ModelPainter painter;

    try
    {
        painter.AddDataObject(nullptr);
        FAIL() << "Expected AddDataObject(nullptr) to throw.";
    }
    catch (const std::runtime_error & ex)
    {
        EXPECT_NE(std::string(ex.what()).find("AddDataObject"), std::string::npos);
    }

    try
    {
        painter.AddReferenceDataObject(nullptr, "ref");
        FAIL() << "Expected AddReferenceDataObject(nullptr, ...) to throw.";
    }
    catch (const std::runtime_error & ex)
    {
        EXPECT_NE(std::string(ex.what()).find("AddReferenceDataObject"), std::string::npos);
    }
}

TEST(DataObjectDispatchIterationArchitectureTest, NormalizeMapObjectNormalizesMapValues)
{
    auto map{ MakeMapObject() };
    const auto original_value{ map.GetMapValue(0) };
    const auto original_sd{ map.GetMapValueSD() };
    ASSERT_GT(original_sd, 0.0f);

    rg::NormalizeMapObject(map);

    EXPECT_NEAR(map.GetMapValue(0), original_value / original_sd, 1.0e-5f);
}

TEST(DataObjectDispatchIterationArchitectureTest, PrepareModelObjectSelectsAndInitializesLocalEntries)
{
    auto model{ MakeModelWithBond() };
    for (auto & atom : model->GetAtomList()) atom->SetSelectedFlag(false);
    for (auto & bond : model->GetBondList()) bond->SetSelectedFlag(false);
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
    for (const auto * atom : model->GetSelectedAtomList())
    {
        ASSERT_NE(atom, nullptr);
        EXPECT_NE(atom->GetLocalPotentialEntry(), nullptr);
    }
    for (const auto * bond : model->GetSelectedBondList())
    {
        ASSERT_NE(bond, nullptr);
        EXPECT_NE(bond->GetLocalPotentialEntry(), nullptr);
    }
}

TEST(DataObjectDispatchIterationArchitectureTest, PrepareModelObjectSelectsAllWhenRequested)
{
    auto model{ MakeModelWithBond() };
    rg::ModelPreparationOptions options;
    options.select_all_atoms = true;
    options.select_all_bonds = true;
    options.update_model = true;

    EXPECT_NO_THROW(rg::PrepareModelObject(*model, options));
    EXPECT_EQ(model->GetNumberOfSelectedAtom(), model->GetNumberOfAtom());
    EXPECT_EQ(model->GetNumberOfSelectedBond(), model->GetNumberOfBond());
}

TEST(DataObjectDispatchIterationArchitectureTest, ApplyModelSelectionSelectsByAtomSelectorRules)
{
    auto model{ MakeModelWithBond() };
    auto & atom_list{ model->GetAtomList() };
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

TEST(DataObjectDispatchIterationArchitectureTest, DataObjectDispatchAsHelpersResolveModelAndMap)
{
    auto model{ MakeModelWithBond() };
    auto map{ MakeMapObject() };
    rg::AtomObject atom;

    EXPECT_EQ(rg::AsModelObject(*model), model.get());
    EXPECT_EQ(rg::AsMapObject(*model), nullptr);
    EXPECT_EQ(rg::AsMapObject(map), &map);
    EXPECT_EQ(rg::AsModelObject(map), nullptr);
    EXPECT_EQ(rg::AsModelObject(atom), nullptr);
    EXPECT_EQ(rg::AsMapObject(atom), nullptr);

    const rg::DataObjectBase & const_model_ref{ *model };
    const rg::DataObjectBase & const_map_ref{ map };
    EXPECT_EQ(rg::AsModelObject(const_model_ref), model.get());
    EXPECT_EQ(rg::AsMapObject(const_map_ref), &map);
}

TEST(DataObjectDispatchIterationArchitectureTest, DataObjectDispatchExpectHelpersResolveModelAndMap)
{
    auto model{ MakeModelWithBond() };
    auto map{ MakeMapObject() };

    const auto & model_ref{ rg::ExpectModelObject(*model, "dispatch-test-model") };
    const auto & map_ref{ rg::ExpectMapObject(map, "dispatch-test-map") };
    EXPECT_EQ(&model_ref, model.get());
    EXPECT_EQ(&map_ref, &map);
}

TEST(DataObjectDispatchIterationArchitectureTest, DataObjectDispatchThrowsOnUnsupportedTargetType)
{
    rg::AtomObject atom;
    EXPECT_THROW(
        (void)rg::ExpectModelObject(atom, "dispatch-test-model"),
        std::runtime_error);
    EXPECT_THROW(
        (void)rg::ExpectMapObject(atom, "dispatch-test-map"),
        std::runtime_error);
}

TEST(DataObjectDispatchIterationArchitectureTest, DataObjectDispatchCatalogTypeNameUsesStableTopLevelNames)
{
    auto model{ MakeModelWithBond() };
    auto map{ MakeMapObject() };
    rg::AtomObject atom;

    EXPECT_EQ(rg::GetCatalogTypeName(*model), "model");
    EXPECT_EQ(rg::GetCatalogTypeName(map), "map");
    EXPECT_THROW((void)rg::GetCatalogTypeName(atom), std::runtime_error);
}

TEST(DataObjectDispatchIterationArchitectureTest, FileProcessFactoryOutputRejectsWrongObjectTypeThroughDispatch)
{
    rg::ModelObjectFactory model_factory;
    rg::MapObjectFactory map_factory;
    rg::AtomObject atom;

    EXPECT_THROW(
        model_factory.OutputDataObject("unused.pdb", atom),
        std::runtime_error);
    EXPECT_THROW(
        map_factory.OutputDataObject("unused.map", atom),
        std::runtime_error);
}

TEST(DataObjectDispatchIterationArchitectureTest, CollectModelAtomsSupportsSelectionAndEntryFilters)
{
    auto model{ MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    ASSERT_EQ(atoms.size(), 2);
    atoms[0]->SetSelectedFlag(true);
    atoms[1]->SetSelectedFlag(false);
    atoms[0]->AddLocalPotentialEntry(std::make_unique<rg::LocalPotentialEntry>());
    model->Update();

    rg::ModelAtomCollectorOptions selected_only_options;
    selected_only_options.selected_only = true;
    selected_only_options.require_local_potential_entry = false;
    const auto selected_only_atoms{ rg::CollectModelAtoms(*model, selected_only_options) };
    ASSERT_EQ(selected_only_atoms.size(), 1);
    EXPECT_EQ(selected_only_atoms.front(), atoms[0].get());

    rg::ModelAtomCollectorOptions require_entry_options;
    require_entry_options.selected_only = false;
    require_entry_options.require_local_potential_entry = true;
    const auto require_entry_atoms{ rg::CollectModelAtoms(*model, require_entry_options) };
    ASSERT_EQ(require_entry_atoms.size(), 1);
    EXPECT_EQ(require_entry_atoms.front(), atoms[0].get());
}

TEST(DataObjectDispatchIterationArchitectureTest, PrepareSimulationAtomsCollectsAtomChargeAndRange)
{
    auto model{ MakeModelWithBond() };
    rg::SimulationAtomPreparationOptions options;
    options.partial_charge_choice = rg::PartialCharge::NEUTRAL;
    options.include_unknown_atoms = true;

    const auto result{ rg::PrepareSimulationAtoms(*model, options) };
    ASSERT_TRUE(result.has_atom);
    EXPECT_EQ(result.atom_list.size(), 2);
    ASSERT_EQ(result.atom_charge_map.size(), 2);
    EXPECT_DOUBLE_EQ(result.atom_charge_map.at(1), 0.0);
    EXPECT_DOUBLE_EQ(result.atom_charge_map.at(2), 0.0);

    const auto range_min{ result.range_minimum };
    const auto range_max{ result.range_maximum };
    EXPECT_FLOAT_EQ(range_min[0], 0.0f);
    EXPECT_FLOAT_EQ(range_max[0], 1.0f);
}

TEST(DataObjectDispatchIterationArchitectureTest, BuildModelAtomBondContextBuildsSelectedContextMaps)
{
    auto model{ MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    auto & bonds{ model->GetBondList() };
    ASSERT_EQ(atoms.size(), 2);
    ASSERT_EQ(bonds.size(), 1);
    atoms[0]->SetSelectedFlag(true);
    atoms[1]->SetSelectedFlag(true);
    bonds[0]->SetSelectedFlag(true);
    model->Update();

    const auto context{ rg::BuildModelAtomBondContext(*model) };
    ASSERT_EQ(context.atom_map.size(), 2);
    ASSERT_EQ(context.bond_map.size(), 2);
    EXPECT_EQ(context.atom_map.at(1), atoms[0].get());
    EXPECT_EQ(context.atom_map.at(2), atoms[1].get());
    EXPECT_EQ(context.bond_map.at(1).size(), 1);
    EXPECT_EQ(context.bond_map.at(2).size(), 1);
}
