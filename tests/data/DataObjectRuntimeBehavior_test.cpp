#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/ModelDerivedState.hpp"
#include "data/detail/ModelObjectParts.hpp"
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
        rg::ModelAnalysisData::Of(*model).Atoms().EnsureLocalEntry(*atom);
    }
    for (auto * bond : model->GetSelectedBonds())
    {
        rg::ModelAnalysisData::Of(*model).Bonds().EnsureLocalEntry(*bond);
    }

    EXPECT_EQ(model->GetSelectedAtomCount(), model->GetNumberOfAtom());
    EXPECT_EQ(model->GetSelectedBondCount(), model->GetNumberOfBond());
    for (const auto * atom : model->GetSelectedAtoms())
    {
        ASSERT_NE(atom, nullptr);
        EXPECT_NE(rg::ModelAnalysisData::FindLocalEntry(*atom), nullptr);
    }
    for (const auto * bond : model->GetSelectedBonds())
    {
        ASSERT_NE(bond, nullptr);
        EXPECT_NE(rg::ModelAnalysisData::FindLocalEntry(*bond), nullptr);
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
    rg::ModelAnalysisData::Of(*model).Atoms().EnsureLocalEntry(*atoms[0]);

    const auto & selected_only_atoms{ model->GetSelectedAtoms() };
    ASSERT_EQ(selected_only_atoms.size(), 1);
    EXPECT_EQ(selected_only_atoms.front(), atoms[0].get());

    std::vector<rg::AtomObject *> require_entry_atoms;
    for (auto & atom : model->GetAtomList())
    {
        if (rg::ModelAnalysisData::FindLocalEntry(*atom) != nullptr)
        {
            require_entry_atoms.emplace_back(atom.get());
        }
    }
    ASSERT_EQ(require_entry_atoms.size(), 1);
    EXPECT_EQ(require_entry_atoms.front(), atoms[0].get());
}

TEST(DataObjectRuntimeBehaviorTest, ModelAnalysisDataCanClearTransientFitStatesWithoutDroppingEntries)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto * bond{ model->GetBondList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };

    analysis_data.Atoms().EnsureLocalEntry(*atom);
    analysis_data.Bonds().EnsureLocalEntry(*bond);
    analysis_data.Atoms().EnsureGroupEntry("atom_class");
    analysis_data.Bonds().EnsureGroupEntry("bond_class");
    analysis_data.Atoms().EnsureFitState(*atom);
    analysis_data.Bonds().EnsureFitState(*bond);

    ASSERT_NE(analysis_data.Atoms().FindLocalEntry(*atom), nullptr);
    ASSERT_NE(analysis_data.Bonds().FindLocalEntry(*bond), nullptr);
    ASSERT_NE(analysis_data.Atoms().FindGroupEntry("atom_class"), nullptr);
    ASSERT_NE(analysis_data.Bonds().FindGroupEntry("bond_class"), nullptr);
    ASSERT_NE(analysis_data.Atoms().FindFitState(*atom), nullptr);
    ASSERT_NE(analysis_data.Bonds().FindFitState(*bond), nullptr);

    analysis_data.ClearTransientFitStates();

    EXPECT_NE(analysis_data.Atoms().FindLocalEntry(*atom), nullptr);
    EXPECT_NE(analysis_data.Bonds().FindLocalEntry(*bond), nullptr);
    EXPECT_NE(analysis_data.Atoms().FindGroupEntry("atom_class"), nullptr);
    EXPECT_NE(analysis_data.Bonds().FindGroupEntry("bond_class"), nullptr);
    EXPECT_EQ(analysis_data.Atoms().FindFitState(*atom), nullptr);
    EXPECT_EQ(analysis_data.Bonds().FindFitState(*bond), nullptr);
}

TEST(DataObjectRuntimeBehaviorTest, ModelAnalysisDataClearDropsEntriesAndFitStates)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto * bond{ model->GetBondList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };

    analysis_data.Atoms().EnsureLocalEntry(*atom);
    analysis_data.Bonds().EnsureLocalEntry(*bond);
    analysis_data.Atoms().EnsureGroupEntry("atom_class");
    analysis_data.Bonds().EnsureGroupEntry("bond_class");
    analysis_data.Atoms().EnsureFitState(*atom);
    analysis_data.Bonds().EnsureFitState(*bond);

    analysis_data.Clear();

    EXPECT_EQ(analysis_data.Atoms().FindLocalEntry(*atom), nullptr);
    EXPECT_EQ(analysis_data.Bonds().FindLocalEntry(*bond), nullptr);
    EXPECT_EQ(analysis_data.Atoms().FindGroupEntry("atom_class"), nullptr);
    EXPECT_EQ(analysis_data.Bonds().FindGroupEntry("bond_class"), nullptr);
    EXPECT_EQ(analysis_data.Atoms().FindFitState(*atom), nullptr);
    EXPECT_EQ(analysis_data.Bonds().FindFitState(*bond), nullptr);
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

TEST(DataObjectRuntimeBehaviorTest, AssemblyBuildsSerialIndexWithoutUsingMovedFromPointer)
{
    static_assert(std::is_default_constructible_v<rg::ModelObjectParts>);
    static_assert(!std::is_copy_constructible_v<rg::ModelObjectParts>);
    static_assert(!std::is_copy_assignable_v<rg::ModelObjectParts>);
    static_assert(std::is_move_constructible_v<rg::ModelObjectParts>);
    static_assert(std::is_move_assignable_v<rg::ModelObjectParts>);

    rg::ModelObjectParts parts;
    auto atom{ std::make_unique<rg::AtomObject>() };
    atom->SetSerialID(42);
    atom->SetPosition(3.0f, 4.0f, 5.0f);
    auto * atom_ptr{ atom.get() };

    parts.atom_list.emplace_back(std::move(atom));
    auto model{ rg::AssembleModelObject(std::move(parts)) };
    EXPECT_EQ(model.GetNumberOfAtom(), 1);
    EXPECT_EQ(model.GetAtomList().front().get(), atom_ptr);
    ASSERT_NE(model.FindAtomPtr(42), nullptr);
    EXPECT_FLOAT_EQ(model.GetCenterOfMassPosition().at(0), 3.0f);
}

TEST(DataObjectRuntimeBehaviorTest, AssemblyBuildsModelWithFreshIndicesAndDerivedCaches)
{
    rg::ModelObjectParts parts;
    parts.atom_list.reserve(2);
    auto atom_1{ std::make_unique<rg::AtomObject>() };
    atom_1->SetSerialID(11);
    atom_1->SetPosition(5.0f, 0.0f, 0.0f);
    auto atom_2{ std::make_unique<rg::AtomObject>() };
    atom_2->SetSerialID(12);
    atom_2->SetPosition(9.0f, 0.0f, 0.0f);
    parts.atom_list.emplace_back(std::move(atom_1));
    parts.atom_list.emplace_back(std::move(atom_2));

    auto model{
        std::make_unique<rg::ModelObject>(rg::AssembleModelObject(std::move(parts))) };
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

TEST(DataObjectRuntimeBehaviorTest, AssemblyInitializesOwnersSelectionAndDerivedState)
{
    rg::ModelObjectParts parts;
    parts.atom_list.reserve(2);

    auto atom_1{ std::make_unique<rg::AtomObject>() };
    atom_1->SetSerialID(1);
    atom_1->SetPosition(0.0f, 0.0f, 0.0f);
    auto atom_2{ std::make_unique<rg::AtomObject>() };
    atom_2->SetSerialID(2);
    atom_2->SetPosition(1.0f, 0.0f, 0.0f);

    auto * atom_1_ptr{ atom_1.get() };
    auto * atom_2_ptr{ atom_2.get() };
    parts.atom_list.emplace_back(std::move(atom_1));
    parts.atom_list.emplace_back(std::move(atom_2));
    auto bond{ std::make_unique<rg::BondObject>(atom_1_ptr, atom_2_ptr) };
    auto * bond_ptr{ bond.get() };
    parts.bond_list.emplace_back(std::move(bond));

    auto model{ rg::AssembleModelObject(std::move(parts)) };

    ASSERT_EQ(model.GetSelectedAtomCount(), 0);
    ASSERT_EQ(model.GetSelectedBondCount(), 0);
    ASSERT_EQ(model.GetNumberOfBond(), 1);
    EXPECT_EQ(model.GetBondList().at(0).get(), bond_ptr);
    EXPECT_EQ(rg::ModelAnalysisData::OwnerOf(*model.GetAtomList().at(0)), &model);
    EXPECT_EQ(rg::ModelAnalysisData::OwnerOf(*model.GetBondList().at(0)), &model);
    EXPECT_FLOAT_EQ(model.GetCenterOfMassPosition().at(0), 0.5f);
    const auto x_range{ model.GetModelPositionRange(0) };
    EXPECT_DOUBLE_EQ(std::get<0>(x_range), 0.0);
    EXPECT_DOUBLE_EQ(std::get<1>(x_range), 1.0);
}

TEST(DataObjectRuntimeBehaviorTest, DerivedStateCachesGeometryQueriesAcrossRepeatedReads)
{
    auto model{ data_test::MakeModelWithBond() };
    const auto & atoms{ model->GetAtomList() };
    ASSERT_EQ(atoms.size(), 2);

    const auto center_of_mass_1{ model->GetCenterOfMassPosition() };
    const auto center_of_mass_2{ model->GetCenterOfMassPosition() };
    EXPECT_EQ(center_of_mass_1, center_of_mass_2);

    const auto x_range_1{ model->GetModelPositionRange(0) };
    const auto x_range_2{ model->GetModelPositionRange(0) };
    EXPECT_EQ(x_range_1, x_range_2);

    const auto first_hits{
        rg::ModelDerivedState::Of(*model).FindAtomsInRange(*model, *atoms.at(0), 0.25) };
    ASSERT_EQ(first_hits.size(), 1);
    EXPECT_EQ(first_hits.front(), atoms.at(0).get());

    const auto repeated_hits{
        rg::ModelDerivedState::Of(*model).FindAtomsInRange(*model, *atoms.at(0), 1.5) };
    const auto repeated_hits_again{
        rg::ModelDerivedState::Of(*model).FindAtomsInRange(*model, *atoms.at(0), 1.5) };
    ASSERT_EQ(repeated_hits.size(), 2);
    ASSERT_EQ(repeated_hits_again.size(), 2);
    std::vector<int> repeated_serials;
    repeated_serials.reserve(repeated_hits.size());
    for (auto * atom : repeated_hits)
    {
        repeated_serials.emplace_back(atom->GetSerialID());
    }
    std::sort(repeated_serials.begin(), repeated_serials.end());

    std::vector<int> repeated_serials_again;
    repeated_serials_again.reserve(repeated_hits_again.size());
    for (auto * atom : repeated_hits_again)
    {
        repeated_serials_again.emplace_back(atom->GetSerialID());
    }
    std::sort(repeated_serials_again.begin(), repeated_serials_again.end());
    EXPECT_EQ(repeated_serials, repeated_serials_again);
}

TEST(DataObjectRuntimeBehaviorTest, DerivedStateSpatialQueriesRemainAvailableAfterAssemblyCopyAndMove)
{
    rg::ModelObjectParts parts;
    parts.atom_list.reserve(2);
    auto atom_1{ std::make_unique<rg::AtomObject>() };
    atom_1->SetSerialID(41);
    atom_1->SetPosition(0.0f, 0.0f, 0.0f);
    auto atom_2{ std::make_unique<rg::AtomObject>() };
    atom_2->SetSerialID(42);
    atom_2->SetPosition(1.0f, 0.0f, 0.0f);
    parts.atom_list.emplace_back(std::move(atom_1));
    parts.atom_list.emplace_back(std::move(atom_2));

    auto assembled_model{ rg::AssembleModelObject(std::move(parts)) };
    const auto assembled_hits{
        rg::ModelDerivedState::Of(assembled_model).FindAtomsInRange(
            assembled_model,
            *assembled_model.GetAtomList().at(0),
            1.5) };
    ASSERT_EQ(assembled_hits.size(), 2);

    auto copied_model{ assembled_model };
    const auto copied_hits{
        rg::ModelDerivedState::Of(copied_model).FindAtomsInRange(
            copied_model,
            *copied_model.GetAtomList().at(0),
            1.5) };
    ASSERT_EQ(copied_hits.size(), 2);
    std::vector<int> copied_serials;
    copied_serials.reserve(copied_hits.size());
    for (auto * atom : copied_hits)
    {
        copied_serials.emplace_back(atom->GetSerialID());
    }
    std::sort(copied_serials.begin(), copied_serials.end());
    EXPECT_EQ(copied_serials, (std::vector<int>{ 41, 42 }));

    auto moved_model{ std::move(copied_model) };
    const auto moved_hits{
        rg::ModelDerivedState::Of(moved_model).FindAtomsInRange(
            moved_model,
            *moved_model.GetAtomList().at(0),
            1.5) };
    ASSERT_EQ(moved_hits.size(), 2);
    std::vector<int> moved_serials;
    moved_serials.reserve(moved_hits.size());
    for (auto * atom : moved_hits)
    {
        moved_serials.emplace_back(atom->GetSerialID());
    }
    std::sort(moved_serials.begin(), moved_serials.end());
    EXPECT_EQ(moved_serials, (std::vector<int>{ 41, 42 }));
}
