#include <gtest/gtest.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/ModelDerivedState.hpp"
#include "data/detail/ModelObjectParts.hpp"
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

namespace {

std::unique_ptr<rg::ModelObject> MakeSpatialQueryModel()
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.reserve(4);

    auto atom_1{ std::make_unique<rg::AtomObject>() };
    atom_1->SetSerialID(1);
    atom_1->SetPosition(0.0f, 0.0f, 0.0f);

    auto atom_2{ std::make_unique<rg::AtomObject>() };
    atom_2->SetSerialID(2);
    atom_2->SetPosition(1.0f, 0.0f, 0.0f);

    auto atom_3{ std::make_unique<rg::AtomObject>() };
    atom_3->SetSerialID(3);
    atom_3->SetPosition(-1.0f, 0.0f, 0.0f);

    auto atom_4{ std::make_unique<rg::AtomObject>() };
    atom_4->SetSerialID(4);
    atom_4->SetPosition(3.0f, 0.0f, 0.0f);

    atom_list.emplace_back(std::move(atom_1));
    atom_list.emplace_back(std::move(atom_2));
    atom_list.emplace_back(std::move(atom_3));
    atom_list.emplace_back(std::move(atom_4));
    return std::make_unique<rg::ModelObject>(std::move(atom_list));
}

std::vector<int> CollectSerialIds(const std::vector<rg::AtomObject *> & atom_list)
{
    std::vector<int> serial_ids;
    serial_ids.reserve(atom_list.size());
    for (const auto * atom : atom_list)
    {
        serial_ids.emplace_back(atom->GetSerialID());
    }
    return serial_ids;
}

} // namespace

TEST(DataObjectAssemblySpatialQueryTest, AssemblyBuildsSerialIndexWithoutUsingMovedFromPointer)
{
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

TEST(DataObjectAssemblySpatialQueryTest, AssemblyBuildsModelWithFreshIndicesAndDerivedCaches)
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

TEST(DataObjectAssemblySpatialQueryTest, AssemblyInitializesOwnersSelectionAndDerivedState)
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

TEST(DataObjectAssemblySpatialQueryTest, DerivedStateCachesGeometryQueriesAcrossRepeatedReads)
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

TEST(DataObjectAssemblySpatialQueryTest, DerivedStateSpatialQueriesRemainAvailableAfterAssemblyCopyAndMove)
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

TEST(DataObjectAssemblySpatialQueryTest, PublicNeighborQueryReturnsDeterministicDistanceSortedNeighbors)
{
    auto model{ MakeSpatialQueryModel() };
    const auto & atoms{ model->GetAtomList() };
    ASSERT_EQ(atoms.size(), 4);

    const auto neighbor_atoms{ model->FindNeighborAtoms(*atoms.at(0), 1.1) };
    EXPECT_EQ(CollectSerialIds(neighbor_atoms), (std::vector<int>{ 2, 3 }));

    const auto neighbor_atoms_from_atom{ atoms.at(0)->FindNeighborAtoms(1.1) };
    EXPECT_EQ(CollectSerialIds(neighbor_atoms_from_atom), (std::vector<int>{ 2, 3 }));
}

TEST(DataObjectAssemblySpatialQueryTest, PublicNeighborQueryCanOptionallyIncludeCenterAtom)
{
    auto model{ MakeSpatialQueryModel() };
    const auto & atoms{ model->GetAtomList() };

    const auto neighbors_without_center{ model->FindNeighborAtoms(*atoms.at(0), 1.1) };
    EXPECT_EQ(CollectSerialIds(neighbors_without_center), (std::vector<int>{ 2, 3 }));

    const auto neighbors_with_center{ model->FindNeighborAtoms(*atoms.at(0), 1.1, true) };
    EXPECT_EQ(CollectSerialIds(neighbors_with_center), (std::vector<int>{ 1, 2, 3 }));

    const auto neighbors_with_self_from_atom{ atoms.at(0)->FindNeighborAtoms(1.1, true) };
    EXPECT_EQ(CollectSerialIds(neighbors_with_self_from_atom), (std::vector<int>{ 1, 2, 3 }));
}

TEST(DataObjectAssemblySpatialQueryTest, PublicNeighborQueryRejectsInvalidInputs)
{
    auto model{ MakeSpatialQueryModel() };
    auto other_model{ data_test::MakeModelWithBond() };
    const auto & atoms{ model->GetAtomList() };
    const auto & other_atoms{ other_model->GetAtomList() };

    EXPECT_THROW(model->FindNeighborAtoms(*atoms.at(0), -0.1), std::invalid_argument);
    EXPECT_THROW(model->FindNeighborAtoms(*other_atoms.at(0), 1.0), std::invalid_argument);

    auto detached_atom{ *atoms.at(0) };
    EXPECT_THROW(detached_atom.FindNeighborAtoms(1.0), std::runtime_error);
}

TEST(DataObjectAssemblySpatialQueryTest, PublicNeighborQueryRefreshesAfterAtomCoordinateChanges)
{
    auto model{ MakeSpatialQueryModel() };
    const auto & atoms{ model->GetAtomList() };

    const auto initial_neighbors{ model->FindNeighborAtoms(*atoms.at(0), 1.1) };
    EXPECT_EQ(CollectSerialIds(initial_neighbors), (std::vector<int>{ 2, 3 }));

    atoms.at(3)->SetPosition(0.5f, 0.0f, 0.0f);

    const auto refreshed_neighbors{ model->FindNeighborAtoms(*atoms.at(0), 1.1) };
    EXPECT_EQ(CollectSerialIds(refreshed_neighbors), (std::vector<int>{ 4, 2, 3 }));
}
