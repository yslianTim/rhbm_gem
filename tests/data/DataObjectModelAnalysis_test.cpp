#include <gtest/gtest.h>

#include <array>
#include <limits>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(DataObjectModelAnalysisTest, SelectedModelEntriesCanBeInitializedForTypedWorkflows)
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
        rg::ModelAnalysisData::Of(*model).EnsureAtomLocalEntry(*atom);
    }
    for (auto * bond : model->GetSelectedBonds())
    {
        rg::ModelAnalysisData::Of(*model).EnsureBondLocalEntry(*bond);
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

TEST(DataObjectModelAnalysisTest, AtomSelectorCanDriveModelSelectionState)
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

TEST(DataObjectModelAnalysisTest, ModelSelectionAndLocalEntriesRemainDirectlyQueryable)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    ASSERT_EQ(atoms.size(), 2);
    model->SetAtomSelected(1, true);
    model->SetAtomSelected(2, false);
    rg::ModelAnalysisData::Of(*model).EnsureAtomLocalEntry(*atoms[0]);

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

TEST(DataObjectModelAnalysisTest, ModelAnalysisDataCanClearTransientFitStatesWithoutDroppingEntries)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto * bond{ model->GetBondList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };

    auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom) };
    auto & bond_entry{ analysis_data.EnsureBondLocalEntry(*bond) };
    analysis_data.EnsureAtomGroupEntry("atom_class");
    analysis_data.EnsureBondGroupEntry("bond_class");
    atom_entry.SetDataset({});
    atom_entry.SetFitResult({});
    bond_entry.SetDataset({});
    bond_entry.SetFitResult({});

    ASSERT_NE(analysis_data.FindAtomLocalEntry(*atom), nullptr);
    ASSERT_NE(analysis_data.FindBondLocalEntry(*bond), nullptr);
    ASSERT_NE(analysis_data.FindAtomGroupEntry("atom_class"), nullptr);
    ASSERT_NE(analysis_data.FindBondGroupEntry("bond_class"), nullptr);
    ASSERT_TRUE(atom_entry.HasDataset());
    ASSERT_TRUE(atom_entry.HasFitResult());
    ASSERT_TRUE(bond_entry.HasDataset());
    ASSERT_TRUE(bond_entry.HasFitResult());

    analysis_data.ClearTransientFitStates();

    const auto * cleared_atom_entry{ analysis_data.FindAtomLocalEntry(*atom) };
    const auto * cleared_bond_entry{ analysis_data.FindBondLocalEntry(*bond) };
    ASSERT_NE(cleared_atom_entry, nullptr);
    ASSERT_NE(cleared_bond_entry, nullptr);
    EXPECT_NE(analysis_data.FindAtomGroupEntry("atom_class"), nullptr);
    EXPECT_NE(analysis_data.FindBondGroupEntry("bond_class"), nullptr);
    EXPECT_FALSE(cleared_atom_entry->HasDataset());
    EXPECT_FALSE(cleared_atom_entry->HasFitResult());
    EXPECT_FALSE(cleared_bond_entry->HasDataset());
    EXPECT_FALSE(cleared_bond_entry->HasFitResult());
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisDataClearDropsEntriesAndFitStates)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto * bond{ model->GetBondList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };

    auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom) };
    auto & bond_entry{ analysis_data.EnsureBondLocalEntry(*bond) };
    analysis_data.EnsureAtomGroupEntry("atom_class");
    analysis_data.EnsureBondGroupEntry("bond_class");
    atom_entry.SetDataset({});
    atom_entry.SetFitResult({});
    bond_entry.SetDataset({});
    bond_entry.SetFitResult({});

    analysis_data.Clear();

    EXPECT_EQ(analysis_data.FindAtomLocalEntry(*atom), nullptr);
    EXPECT_EQ(analysis_data.FindBondLocalEntry(*bond), nullptr);
    EXPECT_EQ(analysis_data.FindAtomGroupEntry("atom_class"), nullptr);
    EXPECT_EQ(analysis_data.FindBondGroupEntry("bond_class"), nullptr);
}

TEST(DataObjectModelAnalysisTest, LocalPotentialEntryCanClearTransientFitPayload)
{
    rg::LocalPotentialEntry entry;
    entry.SetDataset({});
    entry.SetFitResult({});

    ASSERT_TRUE(entry.HasDataset());
    ASSERT_TRUE(entry.HasFitResult());

    entry.ClearTransientFitState();

    EXPECT_FALSE(entry.HasDataset());
    EXPECT_FALSE(entry.HasFitResult());
}

TEST(DataObjectModelAnalysisTest, LocalPotentialEntryPreservesFitResultStatus)
{
    rg::LocalPotentialEntry entry;
    HRLBetaEstimateResult fit_result;
    fit_result.status = HRLEstimationStatus::NUMERICAL_FALLBACK;

    entry.SetFitResult(fit_result);

    ASSERT_TRUE(entry.HasFitResult());
    EXPECT_EQ(HRLEstimationStatus::NUMERICAL_FALLBACK, entry.GetFitResult().status);
}

TEST(DataObjectModelAnalysisTest, RebuildAtomGroupEntriesFromSelectionTracksOnlySelectedAtoms)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * first_atom{ model->GetAtomList().at(0).get() };
    auto * second_atom{ model->GetAtomList().at(1).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    const auto & simple_class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };

    model->SelectAllAtoms(false);
    model->SetAtomSelected(first_atom->GetSerialID(), true);
    analysis_data.EnsureAtomGroupEntry("stale_atom_class").AddMember(999, *second_atom);

    analysis_data.RebuildAtomGroupEntriesFromSelection(*model);

    EXPECT_EQ(analysis_data.FindAtomGroupEntry("stale_atom_class"), nullptr);
    const auto * simple_group_entry{ analysis_data.FindAtomGroupEntry(simple_class_key) };
    ASSERT_NE(simple_group_entry, nullptr);

    size_t member_count{ 0 };
    for (const auto group_key : simple_group_entry->CollectGroupKeys())
    {
        member_count += simple_group_entry->GetMemberCount(group_key);
        for (const auto * atom : simple_group_entry->GetMembers(group_key))
        {
            ASSERT_NE(atom, nullptr);
            EXPECT_EQ(atom, first_atom);
        }
    }
    EXPECT_EQ(member_count, 1U);

    model->SetAtomSelected(first_atom->GetSerialID(), false);
    model->SetAtomSelected(second_atom->GetSerialID(), true);
    analysis_data.RebuildAtomGroupEntriesFromSelection(*model);

    simple_group_entry = analysis_data.FindAtomGroupEntry(simple_class_key);
    ASSERT_NE(simple_group_entry, nullptr);
    member_count = 0;
    for (const auto group_key : simple_group_entry->CollectGroupKeys())
    {
        member_count += simple_group_entry->GetMemberCount(group_key);
        for (const auto * atom : simple_group_entry->GetMembers(group_key))
        {
            ASSERT_NE(atom, nullptr);
            EXPECT_EQ(atom, second_atom);
        }
    }
    EXPECT_EQ(member_count, 1U);
}

TEST(DataObjectModelAnalysisTest, RebuildBondGroupEntriesFromSelectionTracksOnlySelectedBonds)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * bond{ model->GetBondList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    const auto & simple_class_key{ ChemicalDataHelper::GetSimpleBondClassKey() };

    model->SelectAllBonds(false);
    analysis_data.EnsureBondGroupEntry("stale_bond_class").AddMember(888, *bond);
    analysis_data.RebuildBondGroupEntriesFromSelection(*model);

    EXPECT_EQ(analysis_data.FindBondGroupEntry("stale_bond_class"), nullptr);
    const auto * simple_group_entry{ analysis_data.FindBondGroupEntry(simple_class_key) };
    ASSERT_NE(simple_group_entry, nullptr);
    size_t member_count{ 0 };
    for (const auto group_key : simple_group_entry->CollectGroupKeys())
    {
        member_count += simple_group_entry->GetMemberCount(group_key);
    }
    EXPECT_EQ(member_count, 0U);

    model->SetBondSelected(
        bond->GetAtomSerialID1(),
        bond->GetAtomSerialID2(),
        true);
    analysis_data.RebuildBondGroupEntriesFromSelection(*model);

    simple_group_entry = analysis_data.FindBondGroupEntry(simple_class_key);
    ASSERT_NE(simple_group_entry, nullptr);
    member_count = 0;
    for (const auto group_key : simple_group_entry->CollectGroupKeys())
    {
        member_count += simple_group_entry->GetMemberCount(group_key);
        for (const auto * group_bond : simple_group_entry->GetMembers(group_key))
        {
            ASSERT_NE(group_bond, nullptr);
            EXPECT_EQ(group_bond, bond);
        }
    }
    EXPECT_EQ(member_count, 1U);
}

TEST(DataObjectModelAnalysisTest, CollectAtomGroupKeysReturnsRebuiltGroupKeySet)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    const auto & simple_class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };

    model->SelectAllAtoms();
    analysis_data.RebuildAtomGroupEntriesFromSelection(*model);

    const auto group_keys{ analysis_data.CollectAtomGroupKeys(simple_class_key) };
    const auto * group_entry{ analysis_data.FindAtomGroupEntry(simple_class_key) };
    ASSERT_NE(group_entry, nullptr);
    EXPECT_EQ(group_keys.size(), group_entry->GroupCount());
    for (const auto & group_key : group_keys)
    {
        EXPECT_TRUE(group_entry->HasGroup(group_key));
    }

    EXPECT_TRUE(analysis_data.CollectAtomGroupKeys("missing_atom_class").empty());
}

TEST(DataObjectModelAnalysisTest, CollectBondGroupKeysReturnsRebuiltGroupKeySet)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    const auto & simple_class_key{ ChemicalDataHelper::GetSimpleBondClassKey() };

    model->SelectAllBonds();
    analysis_data.RebuildBondGroupEntriesFromSelection(*model);

    const auto group_keys{ analysis_data.CollectBondGroupKeys(simple_class_key) };
    const auto * group_entry{ analysis_data.FindBondGroupEntry(simple_class_key) };
    ASSERT_NE(group_entry, nullptr);
    EXPECT_EQ(group_keys.size(), group_entry->GroupCount());
    for (const auto & group_key : group_keys)
    {
        EXPECT_TRUE(group_entry->HasGroup(group_key));
    }

    EXPECT_TRUE(analysis_data.CollectBondGroupKeys("missing_bond_class").empty());
}

TEST(DataObjectModelAnalysisTest, ModelAtomsExposeStableSerialAndPositionInputsForTypedWorkflows)
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

TEST(DataObjectModelAnalysisTest, SelectedAtomsAndBondsRemainQueryableForContextAssembly)
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
