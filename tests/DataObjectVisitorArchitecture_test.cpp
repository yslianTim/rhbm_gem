#include <gtest/gtest.h>

#include <array>
#include <chrono>
#include <condition_variable>
#include <exception>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "CommandTestHelpers.hpp"
#include "DataObjectManager.hpp"
#include "DataObjectVisitorBase.hpp"
#include "MapInterpolationVisitor.hpp"
#include "MapObject.hpp"
#include "ModelObject.hpp"
#include "SamplerBase.hpp"

namespace rg = rhbm_gem;

namespace {

class EventVisitor : public rg::DataObjectVisitorBase
{
public:
    std::vector<std::string> events;

    void VisitAtomObject(rg::AtomObject * data_object) override
    {
        (void)data_object;
        events.emplace_back("Atom");
    }

    void VisitBondObject(rg::BondObject * data_object) override
    {
        (void)data_object;
        events.emplace_back("Bond");
    }

    void VisitModelObject(rg::ModelObject * data_object) override
    {
        (void)data_object;
        events.emplace_back("Model");
    }

    void VisitMapObject(rg::MapObject * data_object) override
    {
        (void)data_object;
        events.emplace_back("Map");
    }
};

class ModelKeyVisitor : public rg::DataObjectVisitorBase
{
public:
    std::vector<std::string> keys;

    void VisitModelObject(rg::ModelObject * data_object) override
    {
        keys.push_back(data_object->GetKeyTag());
    }
};

class BlockingModelVisitor : public rg::DataObjectVisitorBase
{
public:
    void VisitModelObject(rg::ModelObject * data_object) override
    {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_started = true;
            m_key = data_object->GetKeyTag();
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

class EmptyStrictVisitor : public rg::StrictDataObjectVisitorBase
{
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

TEST(DataObjectVisitorArchitectureTest, ConcreteAcceptDispatchesToMatchingVisitMethod)
{
    auto model{ MakeModelWithBond() };
    auto map{ MakeMapObject() };
    rg::AtomObject atom;
    rg::BondObject bond;
    EventVisitor visitor;

    atom.Accept(&visitor);
    bond.Accept(&visitor);
    map.Accept(&visitor);
    model->Accept(&visitor, rg::ModelVisitMode::SelfOnly);

    EXPECT_EQ(visitor.events, (std::vector<std::string>{ "Atom", "Bond", "Map", "Model" }));
}

TEST(DataObjectVisitorArchitectureTest, ModelLegacyAcceptVisitsAtomsThenModel)
{
    auto model{ MakeModelWithBond() };
    EventVisitor visitor;

    model->Accept(&visitor);

    EXPECT_EQ(visitor.events, (std::vector<std::string>{ "Atom", "Atom", "Model" }));
}

TEST(DataObjectVisitorArchitectureTest, ModelVisitModeBondsThenSelfVisitsBondsAndModelOnly)
{
    auto model{ MakeModelWithBond() };
    EventVisitor visitor;

    model->Accept(&visitor, rg::ModelVisitMode::BondsThenSelf);

    EXPECT_EQ(visitor.events, (std::vector<std::string>{ "Bond", "Model" }));
}

TEST(DataObjectVisitorArchitectureTest, ModelVisitModeAtomsAndBondsThenSelfVisitsAllInOrder)
{
    auto model{ MakeModelWithBond() };
    EventVisitor visitor;

    model->Accept(&visitor, rg::ModelVisitMode::AtomsAndBondsThenSelf);

    EXPECT_EQ(
        visitor.events,
        (std::vector<std::string>{ "Atom", "Atom", "Bond", "Model" }));
}

TEST(DataObjectVisitorArchitectureTest, ModelVisitModeSelfOnlyVisitsOnlyModel)
{
    auto model{ MakeModelWithBond() };
    EventVisitor visitor;

    model->Accept(&visitor, rg::ModelVisitMode::SelfOnly);

    EXPECT_EQ(visitor.events, (std::vector<std::string>{ "Model" }));
}

TEST(DataObjectVisitorArchitectureTest, NullVisitorThrowsInvalidArgumentOnAllAcceptEntrypoints)
{
    auto model{ MakeModelWithBond() };
    auto map{ MakeMapObject() };
    rg::AtomObject atom;
    rg::BondObject bond;
    rg::DataObjectManager manager;
    manager.ProcessFile(command_test::TestDataPath("test_model.cif"), "model");

    auto * null_visitor{ static_cast<rg::DataObjectVisitorBase *>(nullptr) };
    EXPECT_THROW(atom.Accept(null_visitor), std::invalid_argument);
    EXPECT_THROW(bond.Accept(null_visitor), std::invalid_argument);
    EXPECT_THROW(map.Accept(null_visitor), std::invalid_argument);
    EXPECT_THROW(model->Accept(null_visitor), std::invalid_argument);
    EXPECT_THROW(model->Accept(null_visitor, rg::ModelVisitMode::SelfOnly), std::invalid_argument);
    EXPECT_THROW(manager.Accept(null_visitor), std::invalid_argument);
}

TEST(DataObjectVisitorArchitectureTest, ManagerAcceptDeterministicOrderSortsKeysWhenKeyListIsEmpty)
{
    rg::DataObjectManager manager;
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    manager.ProcessFile(model_path, "b_model");
    manager.ProcessFile(model_path, "a_model");

    rg::DataObjectManager::VisitOptions options;
    options.deterministic_order = true;
    options.model_visit_mode = rg::ModelVisitMode::SelfOnly;

    ModelKeyVisitor visitor;
    manager.Accept(&visitor, {}, options);

    EXPECT_EQ(visitor.keys, (std::vector<std::string>{ "a_model", "b_model" }));
}

TEST(DataObjectVisitorArchitectureTest, ManagerAcceptPreservesInputKeyOrderWhenKeyListProvided)
{
    rg::DataObjectManager manager;
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    manager.ProcessFile(model_path, "b_model");
    manager.ProcessFile(model_path, "a_model");

    rg::DataObjectManager::VisitOptions options;
    options.deterministic_order = true;
    options.model_visit_mode = rg::ModelVisitMode::SelfOnly;

    ModelKeyVisitor visitor;
    manager.Accept(&visitor, { "b_model", "a_model" }, options);

    EXPECT_EQ(visitor.keys, (std::vector<std::string>{ "b_model", "a_model" }));
}

TEST(DataObjectVisitorArchitectureTest, ManagerAcceptUsesSnapshotSoClearCanRunConcurrently)
{
    rg::DataObjectManager manager;
    manager.ProcessFile(command_test::TestDataPath("test_model.cif"), "model");

    rg::DataObjectManager::VisitOptions options;
    options.model_visit_mode = rg::ModelVisitMode::SelfOnly;

    BlockingModelVisitor visitor;
    std::exception_ptr worker_error;
    std::thread worker(
        [&]
        {
            try
            {
                manager.Accept(&visitor, { "model" }, options);
            }
            catch (...)
            {
                worker_error = std::current_exception();
            }
        });

    ASSERT_TRUE(visitor.WaitStarted(std::chrono::seconds(2)));
    manager.ClearDataObjects();
    visitor.Release();
    worker.join();

    EXPECT_EQ(worker_error, nullptr);
    EXPECT_EQ(visitor.Key(), "model");
}

TEST(DataObjectVisitorArchitectureTest, StrictVisitorThrowsWhenTypeHandlerIsNotOverridden)
{
    rg::AtomObject atom;
    EmptyStrictVisitor visitor;

    EXPECT_THROW(atom.Accept(&visitor), std::logic_error);
}

TEST(DataObjectVisitorArchitectureTest, MapInterpolationVisitorFailsFastOnUnsupportedDataObjectType)
{
    SinglePointSampler sampler;
    rg::MapInterpolationVisitor visitor{ &sampler };
    rg::AtomObject atom;

    EXPECT_THROW(atom.Accept(&visitor), std::logic_error);
}

TEST(DataObjectVisitorArchitectureTest, MapInterpolationVisitorConsumeSamplingDataClearsState)
{
    auto map{ MakeMapObject() };
    SinglePointSampler sampler;
    rg::MapInterpolationVisitor visitor{ &sampler };
    visitor.SetPosition({ 0.0f, 0.0f, 0.0f });
    visitor.SetAxisVector({ 1.0f, 0.0f, 0.0f });

    map.Accept(&visitor);
    ASSERT_FALSE(visitor.GetSamplingDataList().empty());

    auto sampled_data{ visitor.ConsumeSamplingDataList() };
    EXPECT_FALSE(sampled_data.empty());
    EXPECT_TRUE(visitor.GetSamplingDataList().empty());
}

TEST(DataObjectVisitorArchitectureTest, MapInterpolationVisitorNullMapVisitClearsStaleOutput)
{
    auto map{ MakeMapObject() };
    SinglePointSampler sampler;
    rg::MapInterpolationVisitor visitor{ &sampler };
    visitor.SetPosition({ 0.0f, 0.0f, 0.0f });
    visitor.SetAxisVector({ 1.0f, 0.0f, 0.0f });

    map.Accept(&visitor);
    ASSERT_FALSE(visitor.GetSamplingDataList().empty());

    visitor.VisitMapObject(nullptr);
    EXPECT_TRUE(visitor.GetSamplingDataList().empty());
}
