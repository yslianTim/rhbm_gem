#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "core/detail/LocalPotentialAccess.hpp"
#include "data/detail/ModelObjectBuilder.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "data/detail/MapSpatialIndex.hpp"
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(DataObjectRuntimeBehaviorTest, NormalizeMapObjectNormalizesMapValues)
{
    auto map{ data_test::MakeMapObject() };
    const auto original_value{ map.GetMapValue(0) };
    const auto original_sd{ map.GetMapValueSD() };
    ASSERT_GT(original_sd, 0.0f);

    map.MapValueArrayNormalization();

    EXPECT_NEAR(map.GetMapValue(0), original_value / original_sd, 1.0e-5f);
    EXPECT_NEAR(map.GetMapValueSD(), 1.0f, 1.0e-5f);
}

TEST(DataObjectRuntimeBehaviorTest, SetMapValueArrayRefreshesStatisticsAndPreservesSpatialQueries)
{
    auto map{ data_test::MakeMapObject() };
    std::vector<size_t> initial_hits;
    rg::MapSpatialIndex(map).CollectGridIndicesInRange(
        map.GetGridPosition(0),
        0.1f,
        initial_hits);
    ASSERT_FALSE(initial_hits.empty());

    auto replacement_values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        replacement_values[i] = static_cast<float>(10 + i);
    }

    map.SetMapValueArray(std::move(replacement_values));

    std::vector<size_t> refreshed_hits;
    rg::MapSpatialIndex(map).CollectGridIndicesInRange(
        map.GetGridPosition(0),
        0.1f,
        refreshed_hits);
    EXPECT_EQ(refreshed_hits, initial_hits);
    EXPECT_FLOAT_EQ(map.GetMapValueMin(), 10.0f);
    EXPECT_FLOAT_EQ(map.GetMapValueMax(), 17.0f);
    EXPECT_FLOAT_EQ(map.GetMapValueMean(), 13.5f);
}

TEST(DataObjectRuntimeBehaviorTest, SelectedModelEntriesCanBeInitializedForTypedWorkflows)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms(false);
    model->SelectAllBonds(false);
    ASSERT_EQ(model->GetSelectedAtomCount(), 0);
    ASSERT_EQ(model->GetSelectedBondCount(), 0);

    model->SelectAllAtoms();
    model->SelectAllBonds();
    model->ApplySymmetrySelection(false);
    for (auto * atom : model->GetSelectedAtoms())
    {
        rg::EnsureLocalPotentialEntry(*model, *atom);
    }
    for (auto * bond : model->GetSelectedBonds())
    {
        rg::EnsureLocalPotentialEntry(*model, *bond);
    }

    EXPECT_EQ(model->GetSelectedAtomCount(), model->GetNumberOfAtom());
    EXPECT_EQ(model->GetSelectedBondCount(), model->GetNumberOfBond());
    for (const auto * atom : model->GetSelectedAtoms())
    {
        ASSERT_NE(atom, nullptr);
        EXPECT_NE(rg::FindLocalPotentialEntry(*atom), nullptr);
    }
    for (const auto * bond : model->GetSelectedBonds())
    {
        ASSERT_NE(bond, nullptr);
        EXPECT_NE(rg::FindLocalPotentialEntry(*bond), nullptr);
    }
}

TEST(DataObjectRuntimeBehaviorTest, AtomSelectorCanDriveModelSelectionState)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atom_list{ model->GetAtomList() };
    ASSERT_EQ(atom_list.size(), 2);
    atom_list.at(0)->SetChainID("A");
    atom_list.at(1)->SetChainID("B");

    ::AtomSelector selector;
    selector.PickChainID("A");
    model->SelectAtoms(
        [&selector](const rg::AtomObject & atom)
        {
            return selector.GetSelectionFlag(
                atom.GetChainID(),
                atom.GetResidue(),
                atom.GetElement());
        });

    ASSERT_EQ(model->GetSelectedAtomCount(), 1);
    EXPECT_EQ(model->GetSelectedAtoms().front(), atom_list.at(0).get());
}

TEST(DataObjectRuntimeBehaviorTest, ModelSelectionAndLocalEntriesRemainDirectlyQueryable)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    ASSERT_EQ(atoms.size(), 2);
    model->SetAtomSelected(1, true);
    model->SetAtomSelected(2, false);
    rg::EnsureLocalPotentialEntry(*model, *atoms[0]);

    const auto & selected_only_atoms{ model->GetSelectedAtoms() };
    ASSERT_EQ(selected_only_atoms.size(), 1);
    EXPECT_EQ(selected_only_atoms.front(), atoms[0].get());

    std::vector<rg::AtomObject *> require_entry_atoms;
    for (auto & atom : model->GetAtomList())
    {
        if (rg::FindLocalPotentialEntry(*atom) != nullptr)
        {
            require_entry_atoms.emplace_back(atom.get());
        }
    }
    ASSERT_EQ(require_entry_atoms.size(), 1);
    EXPECT_EQ(require_entry_atoms.front(), atoms[0].get());
}

TEST(DataObjectRuntimeBehaviorTest, ModelAtomsExposeStableSerialAndPositionInputsForTypedWorkflows)
{
    auto model{ data_test::MakeModelWithBond() };
    std::vector<rg::AtomObject *> atom_list;
    std::unordered_map<int, double> atom_charge_map;
    std::array<float, 3> range_minimum{
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max() };
    std::array<float, 3> range_maximum{
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest() };
    for (auto & atom : model->GetAtomList())
    {
        atom_list.emplace_back(atom.get());
        atom_charge_map.emplace(atom->GetSerialID(), 0.0);
        const auto & atom_position{ atom->GetPositionRef() };
        range_minimum[0] = std::min(range_minimum[0], atom_position[0]);
        range_minimum[1] = std::min(range_minimum[1], atom_position[1]);
        range_minimum[2] = std::min(range_minimum[2], atom_position[2]);
        range_maximum[0] = std::max(range_maximum[0], atom_position[0]);
        range_maximum[1] = std::max(range_maximum[1], atom_position[1]);
        range_maximum[2] = std::max(range_maximum[2], atom_position[2]);
    }

    ASSERT_FALSE(atom_list.empty());
    EXPECT_EQ(atom_list.size(), 2);
    ASSERT_EQ(atom_charge_map.size(), 2);
    EXPECT_DOUBLE_EQ(atom_charge_map.at(1), 0.0);
    EXPECT_DOUBLE_EQ(atom_charge_map.at(2), 0.0);
    EXPECT_FLOAT_EQ(range_minimum[0], 0.0f);
    EXPECT_FLOAT_EQ(range_maximum[0], 1.0f);
}

TEST(DataObjectRuntimeBehaviorTest, SelectedAtomsAndBondsRemainQueryableForContextAssembly)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    auto & bonds{ model->GetBondList() };
    ASSERT_EQ(atoms.size(), 2);
    ASSERT_EQ(bonds.size(), 1);
    model->SelectAllAtoms();
    model->SelectAllBonds();

    std::unordered_map<int, rg::AtomObject *> atom_map;
    std::unordered_map<int, std::vector<rg::BondObject *>> bond_map;
    for (auto * atom : model->GetSelectedAtoms())
    {
        atom_map.emplace(atom->GetSerialID(), atom);
    }
    for (auto * bond : model->GetSelectedBonds())
    {
        bond_map[bond->GetAtomSerialID1()].emplace_back(bond);
        bond_map[bond->GetAtomSerialID2()].emplace_back(bond);
    }

    ASSERT_EQ(atom_map.size(), 2);
    ASSERT_EQ(bond_map.size(), 2);
    EXPECT_EQ(atom_map.at(1), atoms[0].get());
    EXPECT_EQ(atom_map.at(2), atoms[1].get());
    EXPECT_EQ(bond_map.at(1).size(), 1);
    EXPECT_EQ(bond_map.at(2).size(), 1);
}

TEST(DataObjectRuntimeBehaviorTest, AddAtomRebuildsSerialIndexWithoutUsingMovedFromPointer)
{
    rg::ModelObjectBuilder builder;
    auto atom{ std::make_unique<rg::AtomObject>() };
    atom->SetSerialID(42);
    atom->SetPosition(3.0f, 4.0f, 5.0f);

    ASSERT_NO_THROW(builder.AddAtom(std::move(atom)));
    auto model{ builder.Build() };
    EXPECT_EQ(model.GetNumberOfAtom(), 1);
    ASSERT_NE(model.FindAtomPtr(42), nullptr);
    EXPECT_FLOAT_EQ(model.GetCenterOfMassPosition().at(0), 3.0f);
}

TEST(DataObjectRuntimeBehaviorTest, BuilderBuildsModelWithFreshIndicesAndDerivedCaches)
{
    rg::ModelObjectBuilder builder;
    std::vector<std::unique_ptr<rg::AtomObject>> replacement_atoms;
    auto atom_1{ std::make_unique<rg::AtomObject>() };
    atom_1->SetSerialID(11);
    atom_1->SetPosition(5.0f, 0.0f, 0.0f);
    auto atom_2{ std::make_unique<rg::AtomObject>() };
    atom_2->SetSerialID(12);
    atom_2->SetPosition(9.0f, 0.0f, 0.0f);
    replacement_atoms.emplace_back(std::move(atom_1));
    replacement_atoms.emplace_back(std::move(atom_2));

    builder.SetAtomList(std::move(replacement_atoms));
    auto model{ std::make_unique<rg::ModelObject>(builder.Build()) };
    model->SetAtomSelected(11, true);
    model->SetAtomSelected(12, false);

    EXPECT_NE(model->FindAtomPtr(11), nullptr);
    EXPECT_NE(model->FindAtomPtr(12), nullptr);
    EXPECT_THROW(model->FindAtomPtr(1), std::out_of_range);
    EXPECT_EQ(model->GetSelectedAtomCount(), 1);
    EXPECT_FLOAT_EQ(model->GetCenterOfMassPosition().at(0), 7.0f);
    const auto x_range{ model->GetModelPositionRange(0) };
    EXPECT_DOUBLE_EQ(std::get<0>(x_range), 5.0);
    EXPECT_DOUBLE_EQ(std::get<1>(x_range), 9.0);
}
