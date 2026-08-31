#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/LocalPotentialSeries.hpp>
#include <rhbm_gem/utils/math/GridSampler.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

namespace rg = rhbm_gem;

namespace {

static_assert(sizeof(double) == 8);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::AtomObject &>().GetPosition()),
    std::array<double, 3>>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::AtomObject &>().GetPositionRef()),
    const std::array<double, 3> &>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::AtomObject &>().GetOccupancy()), double>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::AtomObject &>().GetTemperature()), double>);
static_assert(std::is_same_v<
    decltype(&rg::AtomObject::SetOccupancy),
    void (rg::AtomObject::*)(double)>);
static_assert(std::is_same_v<
    decltype(&rg::AtomObject::SetTemperature),
    void (rg::AtomObject::*)(double)>);
static_assert(std::is_same_v<
    decltype(static_cast<void (rg::AtomObject::*)(
        const std::array<double, 3> &)>(&rg::AtomObject::SetPosition)),
    void (rg::AtomObject::*)(const std::array<double, 3> &)>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::BondObject &>().GetPosition()),
    std::array<double, 3>>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::BondObject &>().GetBondVector()),
    std::array<double, 3>>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::ChemicalComponentEntry &>()
        .GetComponentMolecularWeight()),
    double>);
static_assert(std::is_same_v<
    decltype(&rg::ChemicalComponentEntry::SetComponentMolecularWeight),
    void (rg::ChemicalComponentEntry::*)(double)>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::MapObject &>().GetGridSpacing()),
    std::array<double, 3>>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::MapObject &>().GetOrigin()),
    std::array<double, 3>>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::MapObject &>().GetMapValueArray()),
    const double *>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::MapObject &>().GetMapValue(0)), double>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::MapObject &>().GetMapValueMean()), double>);
static_assert(std::is_same_v<
    decltype(std::declval<const rg::MapObject &>().GetMapValueSD()), double>);
static_assert(std::is_same_v<
    decltype(static_cast<void (rg::MapObject::*)(
        std::unique_ptr<double[]>)>(&rg::MapObject::SetMapValueArray)),
    void (rg::MapObject::*)(std::unique_ptr<double[]>)>);
static_assert(std::is_same_v<
    decltype(std::declval<rg::ModelObject &>().GetCenterOfMassPosition()),
    std::array<double, 3>>);
static_assert(std::is_same_v<decltype(SamplingPoint{}.distance), double>);
static_assert(std::is_same_v<
    decltype(SamplingPoint{}.position), std::array<double, 3>>);
static_assert(std::is_same_v<decltype(LocalPotentialSample{}.response), double>);
static_assert(std::is_same_v<
    decltype(static_cast<void (GridSampler::*)(double)>(
        &GridSampler::SetWindowSize)),
    void (GridSampler::*)(double)>);
static_assert(std::is_same_v<
    decltype(static_cast<void (GridSampler::*)(
        const std::array<double, 3> &)>(
            &GridSampler::SetReferenceUVector)),
    void (GridSampler::*)(const std::array<double, 3> &)>);
static_assert(std::is_same_v<
    decltype(&GridSampler::GenerateSamplingPoints),
    SamplingPointList (GridSampler::*)(
        const std::array<double, 3> &,
        const std::array<double, 3> &) const>);
static_assert(std::is_same_v<
    decltype(&rg::sphere_sampler::GenerateSamplingPointList),
    SamplingPointList (*)(
        const std::array<double, 3> &,
        SphereSamplingMethod)>);
static_assert(std::is_same_v<
    decltype(&rg::sample_filter::FilterSamplingPointList),
    void (*)(
        SamplingPointList &,
        const std::array<double, 3> &,
        const std::vector<std::array<double, 3>> &,
        double)>);
static_assert(std::is_same_v<
    decltype(&rg::core::SampleMapValues),
    LocalPotentialSampleList (*)(
        const rg::MapObject &,
        const GridSampler &,
        const std::array<double, 3> &,
        const std::array<double, 3> &)>);
static_assert(std::is_same_v<
    decltype(rg::local_potential_series::ComputeDistanceRange(
        std::declval<const LocalPotentialSampleList &>())),
    std::tuple<double, double>>);
static_assert(std::is_constructible_v<
    rg::MapObject,
    const std::array<int, 3> &,
    const std::array<double, 3> &,
    const std::array<double, 3> &,
    std::unique_ptr<double[]>>);

TEST(DataObjectPrecisionContractTest, HeldAndDerivedValuesRetainDoublePrecision)
{
    constexpr double beyond_float_precision{ 1.0000000000000002 };
    ASSERT_NE(
        beyond_float_precision,
        static_cast<double>(static_cast<float>(beyond_float_precision)));

    rg::AtomObject atom;
    atom.SetPosition({
        beyond_float_precision,
        16777217.0,
        -0.123456789012345 });
    atom.SetOccupancy(0.123456789012345);
    atom.SetTemperature(12.3456789012345);
    atom.AddAlternatePosition("B", {
        0.987654321098765,
        -16777217.0,
        0.3333333333333333 });
    atom.AddAlternateOccupancy("B", 0.234567890123456);
    atom.AddAlternateTemperature("B", 23.4567890123456);

    const rg::AtomObject atom_copy{ atom };
    EXPECT_DOUBLE_EQ(atom_copy.GetPosition().at(0), beyond_float_precision);
    EXPECT_DOUBLE_EQ(atom_copy.GetPosition().at(1), 16777217.0);
    EXPECT_DOUBLE_EQ(atom_copy.GetOccupancy(), 0.123456789012345);
    EXPECT_DOUBLE_EQ(atom_copy.GetTemperature(), 12.3456789012345);
    EXPECT_DOUBLE_EQ(
        atom_copy.GetAlternatePositions().at("B").at(0),
        0.987654321098765);
    EXPECT_DOUBLE_EQ(
        atom_copy.GetAlternateOccupancies().at("B"),
        0.234567890123456);
    EXPECT_DOUBLE_EQ(
        atom_copy.GetAlternateTemperatures().at("B"),
        23.4567890123456);

    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    atoms.emplace_back(std::make_unique<rg::AtomObject>(atom));
    rg::ModelObject model{ std::move(atoms) };
    EXPECT_EQ(model.GetCenterOfMassPosition(), atom.GetPosition());

    rg::ChemicalComponentEntry component;
    component.SetComponentMolecularWeight(123.456789012345);
    EXPECT_DOUBLE_EQ(
        component.GetComponentMolecularWeight(), 123.456789012345);

    auto map_values{ std::make_unique<double[]>(2) };
    map_values[0] = 16777217.0;
    map_values[1] = 16777218.0;
    rg::MapObject map{
        { 2, 1, 1 },
        { beyond_float_precision, 1.0, 1.0 },
        { 0.123456789012345, 0.0, 0.0 },
        std::move(map_values) };
    const rg::MapObject map_copy{ map };
    EXPECT_DOUBLE_EQ(map_copy.GetGridSpacing().at(0), beyond_float_precision);
    EXPECT_DOUBLE_EQ(map_copy.GetOrigin().at(0), 0.123456789012345);
    EXPECT_DOUBLE_EQ(map_copy.GetMapValue(0), 16777217.0);
    EXPECT_DOUBLE_EQ(map_copy.GetMapValueMean(), 16777217.5);

    const LocalPotentialSampleList samples{
        { 16777217.0, SamplingPoint{ beyond_float_precision } },
        { 16777218.0, SamplingPoint{ 2.0000000000000004 } }
    };
    const auto distance_range{
        rg::local_potential_series::ComputeDistanceRange(samples) };
    const auto response_range{
        rg::local_potential_series::ComputeResponseRange(samples) };
    EXPECT_DOUBLE_EQ(std::get<0>(distance_range), beyond_float_precision);
    EXPECT_DOUBLE_EQ(std::get<1>(distance_range), 2.0000000000000004);
    EXPECT_DOUBLE_EQ(std::get<0>(response_range), 16777217.0);
    EXPECT_DOUBLE_EQ(std::get<1>(response_range), 16777218.0);
}

} // namespace
