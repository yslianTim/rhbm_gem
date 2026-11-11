#include <gtest/gtest.h>
#include "core/sampling/MapInterpolation.hpp"
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace {

struct AtomSpec
{
    std::array<double, 3> position;
    Element element{ Element::CARBON };
};

rhbm_gem::MapObject MakeMapObject(
    const std::array<int, 3> & grid_size,
    const std::array<double, 3> & grid_spacing,
    const std::array<double, 3> & origin,
    const std::vector<double> & values)
{
    auto map_values{ std::make_unique<double[]>(values.size()) };
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        map_values[i] = values[i];
    }
    return rhbm_gem::MapObject{
        grid_size,
        grid_spacing,
        origin,
        std::move(map_values)
    };
}

rhbm_gem::MapObject MakeMapObject(const std::vector<double> & values)
{
    return MakeMapObject(
        {
            static_cast<int>(values.size()),
            1,
            1
        },
        { 1.0, 1.0, 1.0 },
        { 0.0, 0.0, 0.0 },
        values);
}

rhbm_gem::MapObject MakeMapObject(
    const std::array<int, 3> & grid_size,
    const std::array<double, 3> & grid_spacing,
    const std::array<double, 3> & origin,
    const std::function<double(const std::array<double, 3> &)> & value_at_position)
{
    std::vector<double> values;
    values.reserve(
        static_cast<std::size_t>(grid_size.at(0)) *
        static_cast<std::size_t>(grid_size.at(1)) *
        static_cast<std::size_t>(grid_size.at(2)));
    for (int z = 0; z < grid_size.at(2); ++z)
    {
        for (int y = 0; y < grid_size.at(1); ++y)
        {
            for (int x = 0; x < grid_size.at(0); ++x)
            {
                values.emplace_back(
                    value_at_position({
                        origin.at(0) + static_cast<double>(x) * grid_spacing.at(0),
                        origin.at(1) + static_cast<double>(y) * grid_spacing.at(1),
                        origin.at(2) + static_cast<double>(z) * grid_spacing.at(2)
                    }));
            }
        }
    }
    return MakeMapObject(grid_size, grid_spacing, origin, values);
}

rhbm_gem::MapObject MakeGaussianMapObject()
{
    constexpr double sigma{ 0.6 };
    return MakeMapObject(
        { 61, 61, 61 },
        { 0.1, 0.1, 0.1 },
        { -3.0, -3.0, -3.0 },
        [](const std::array<double, 3> & position)
        {
            const auto distance_square{
                position.at(0) * position.at(0) +
                position.at(1) * position.at(1) +
                position.at(2) * position.at(2)
            };
            return std::exp(-0.5 * distance_square / (sigma * sigma));
        });
}

std::unique_ptr<rhbm_gem::ModelObject> MakeModelObject(
    const std::vector<AtomSpec> & atom_specs)
{
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    atoms.reserve(atom_specs.size());
    for (std::size_t i = 0; i < atom_specs.size(); ++i)
    {
        auto atom{ std::make_unique<rhbm_gem::AtomObject>() };
        atom->SetSerialID(static_cast<int>(i + 1));
        atom->SetElement(atom_specs.at(i).element);
        atom->SetPosition(atom_specs.at(i).position);
        atoms.emplace_back(std::move(atom));
    }
    return std::make_unique<rhbm_gem::ModelObject>(std::move(atoms));
}

std::unique_ptr<rhbm_gem::ModelObject> MakeCrowdedModelObject()
{
    std::vector<AtomSpec> atom_specs{
        { { 0.0, 0.0, 0.0 }, Element::CARBON }
    };
    for (int x = -1; x <= 1; ++x)
    {
        for (int y = -1; y <= 1; ++y)
        {
            for (int z = -1; z <= 1; ++z)
            {
                if (x == 0 && y == 0 && z == 0)
                {
                    continue;
                }
                const auto norm{
                    std::sqrt(static_cast<double>(x * x + y * y + z * z))
                };
                atom_specs.emplace_back(AtomSpec{
                    {
                        static_cast<double>(x) / norm,
                        static_cast<double>(y) / norm,
                        static_cast<double>(z) / norm
                    },
                    Element::CARBON
                });
            }
        }
    }
    return MakeModelObject(atom_specs);
}

double ComputeDistanceSquare(
    const std::array<double, 3> & lhs,
    const std::array<double, 3> & rhs)
{
    double distance_square{ 0.0 };
    for (std::size_t axis = 0; axis < lhs.size(); ++axis)
    {
        const auto difference{
            lhs.at(axis) - rhs.at(axis)
        };
        distance_square += difference * difference;
    }
    return distance_square;
}

void ExpectSamplingPointListsEqual(
    const SamplingPointList & lhs,
    const SamplingPointList & rhs)
{
    ASSERT_EQ(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(lhs.at(i).distance, rhs.at(i).distance);
        EXPECT_EQ(lhs.at(i).position, rhs.at(i).position);
        EXPECT_EQ(lhs.at(i).is_selected, rhs.at(i).is_selected);
    }
}

} // namespace

TEST(QScoreHelperTest, ReturnsAmplitudeAndOffsetForTypicalMap)
{
    const auto map{ MakeMapObject({ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 }) };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    constexpr double expected_standard_deviation{ 2.449489742783178 };
    constexpr double expected_offset{ 4.5 - expected_standard_deviation };
    constexpr double expected_height{ 8.0 - expected_offset };
    EXPECT_NEAR(height, expected_height, 1.0e-5);
    EXPECT_NEAR(offset, expected_offset, 1.0e-5);
}

TEST(QScoreHelperTest, CapsReferenceHighAtMeanPlusTenStandardDeviations)
{
    std::vector<double> values(121, 0.0);
    values.back() = 1.0;
    const auto map{ MakeMapObject(values) };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    constexpr double expected_height{ 1.0 / 121.0 + 10.0 / 11.0 };
    EXPECT_NEAR(height, expected_height, 1.0e-5);
    EXPECT_DOUBLE_EQ(offset, 0.0);
    EXPECT_LT(height + offset, map.GetMapValueMax());
}

TEST(QScoreHelperTest, ClampsOffsetToObservedMinimum)
{
    const auto map{ MakeMapObject({ 0.0, 0.0, 0.0, 10.0 }) };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    EXPECT_DOUBLE_EQ(height, 10.0);
    EXPECT_DOUBLE_EQ(offset, 0.0);
}

TEST(QScoreHelperTest, ReturnsZeroAmplitudeForConstantMap)
{
    const auto map{
        MakeMapObject({ 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5 })
    };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    EXPECT_DOUBLE_EQ(height, 0.0);
    EXPECT_DOUBLE_EQ(offset, 3.5);
}

TEST(QScoreHelperTest, ReturnsDeterministicSpiralPointsForIsolatedAtom)
{
    const auto model{
        MakeModelObject({
            { { 1.0, 2.0, 3.0 }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };

    const auto first_points{
        rhbm_gem::core::GetRadialPointsForQScore(atom, *model, 1.0, 4)
    };
    const auto second_points{
        rhbm_gem::core::GetRadialPointsForQScore(atom, *model, 1.0, 4)
    };

    ASSERT_EQ(first_points.size(), 4u);
    ExpectSamplingPointListsEqual(first_points, second_points);
    for (const auto & point : first_points)
    {
        EXPECT_DOUBLE_EQ(point.distance, 1.0);
        EXPECT_TRUE(point.is_selected);
        EXPECT_NEAR(
            ComputeDistanceSquare(point.position, atom.GetPositionRef()),
            1.0,
            1.0e-5);
    }
    EXPECT_NEAR(first_points.front().position.at(0), 1.0, 1.0e-6);
    EXPECT_NEAR(first_points.front().position.at(1), 2.0, 1.0e-6);
    EXPECT_NEAR(first_points.front().position.at(2), 2.0, 1.0e-6);
    EXPECT_NEAR(first_points.back().position.at(0), 1.0, 1.0e-6);
    EXPECT_NEAR(first_points.back().position.at(1), 2.0, 1.0e-6);
    EXPECT_NEAR(first_points.back().position.at(2), 4.0, 1.0e-6);
}

TEST(QScoreHelperTest, ReturnsAllAcceptedPointsFromFirstSuccessfulRetry)
{
    const auto model{
        MakeModelObject({
            { { 0.0, 0.0, 0.0 }, Element::CARBON },
            { { 0.0, 0.0, 0.5 }, Element::CARBON }
        })
    };
    const auto & target_atom{ *model->GetAtomList().at(0) };
    const auto & neighbor_atom{ *model->GetAtomList().at(1) };

    const auto points{
        rhbm_gem::core::GetRadialPointsForQScore(
            target_atom,
            *model,
            1.0,
            2)
    };

    const auto isolated{ MakeModelObject({ { { 0.0, 0.0, 0.0 }, Element::CARBON } }) };
    auto reference{ rhbm_gem::core::GetRadialPointsForQScore(
        *isolated->GetAtomList().front(), *isolated, 1.0, 4) };
    ASSERT_EQ(reference.size(), 4u);
    reference.pop_back(); // The north pole is the only rejected point in the first retry.
    ExpectSamplingPointListsEqual(points, reference);
    ASSERT_EQ(points.size(), 3u);
    for (const auto & point : points)
    {
        EXPECT_GE(
            ComputeDistanceSquare(
                point.position,
                neighbor_atom.GetPositionRef()),
            0.9 * 0.9);
    }
}

TEST(QScoreHelperTest, IgnoresHydrogenAtomsWhenFilteringRadialPoints)
{
    const auto model{
        MakeModelObject({
            { { 0.0, 0.0, 0.0 }, Element::CARBON },
            { { 0.0, 0.0, 0.5 }, Element::HYDROGEN }
        })
    };
    const auto & target_atom{ *model->GetAtomList().at(0) };

    const auto points{
        rhbm_gem::core::GetRadialPointsForQScore(
            target_atom,
            *model,
            1.0,
            2)
    };

    ASSERT_EQ(points.size(), 2u);
    EXPECT_NEAR(points.front().position.at(2), -1.0, 1.0e-6);
    EXPECT_NEAR(points.back().position.at(2), 1.0, 1.0e-6);
}

TEST(QScoreHelperTest, ReturnsEmptyWhenCrowdingPreventsRequestedPointCount)
{
    const auto model{ MakeCrowdedModelObject() };
    const auto & target_atom{ *model->GetAtomList().front() };

    const auto points{
        rhbm_gem::core::GetRadialPointsForQScore(
            target_atom,
            *model,
            1.0,
            8)
    };

    EXPECT_TRUE(points.empty());
}

TEST(QScoreHelperTest, RejectsInvalidRadialPointArguments)
{
    const auto model{
        MakeModelObject({
            { { 0.0, 0.0, 0.0 }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };

    EXPECT_THROW(
        rhbm_gem::core::GetRadialPointsForQScore(atom, *model, 0.0, 2),
        std::invalid_argument);
    EXPECT_THROW(
        rhbm_gem::core::GetRadialPointsForQScore(atom, *model, -1.0, 2),
        std::invalid_argument);
    EXPECT_THROW(
        rhbm_gem::core::GetRadialPointsForQScore(
            atom,
            *model,
            std::numeric_limits<double>::quiet_NaN(),
            2),
        std::invalid_argument);
    EXPECT_THROW(
        rhbm_gem::core::GetRadialPointsForQScore(
            atom,
            *model,
            std::numeric_limits<double>::infinity(),
            2),
        std::invalid_argument);
    EXPECT_THROW(
        rhbm_gem::core::GetRadialPointsForQScore(atom, *model, 1.0, 1),
        std::invalid_argument);
    EXPECT_THROW(
        rhbm_gem::core::GetRadialPointsForQScore(atom, *model, 1.0, -1),
        std::invalid_argument);
}

TEST(QScoreHelperTest, CalculatesHighQScoreForMatchingGaussianWidth)
{
    constexpr double height{ 1.0 };
    constexpr double offset{ 0.0 };
    constexpr double sigma{ 0.4 };
    constexpr double mismatched_sigma{ 0.8 };
    constexpr int point_count{ 8 };
    const auto map{
        MakeMapObject(
            { 41, 41, 41 },
            { 0.1, 0.1, 0.1 },
            { -2.0, -2.0, -2.0 },
            [](const std::array<double, 3> & position)
            {
                const auto distance_square{
                    position.at(0) * position.at(0) +
                    position.at(1) * position.at(1) +
                    position.at(2) * position.at(2)
                };
                return std::exp(-0.5 * distance_square / (sigma * sigma));
            })
    };
    const auto model{
        MakeModelObject({
            { { 0.0, 0.0, 0.0 }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };

    const auto matching_score{
        rhbm_gem::core::CalculateQScoreForAtom(
            atom,
            map,
            *model,
            height,
            offset,
            sigma,
            point_count)
    };
    const auto mismatched_score{
        rhbm_gem::core::CalculateQScoreForAtom(
            atom,
            map,
            *model,
            height,
            offset,
            mismatched_sigma,
            point_count)
    };

    EXPECT_GT(matching_score, 0.99);
    EXPECT_GT(matching_score, mismatched_score);
}

TEST(QScoreHelperTest, UsesCenterWeightMaximumShellAndAllAcceptedPoints)
{
    constexpr double height{ 0.5 };
    constexpr double offset{ 1.0 };
    constexpr double sigma{ 1.0 };
    constexpr int point_count{ 2 };
    const auto map{
        MakeMapObject(
            { 5, 5, 5 },
            { 1.0, 1.0, 1.0 },
            { -2.0, -2.0, -2.0 },
            [](const std::array<double, 3> & position)
            {
                return position.at(2);
            })
    };
    const auto model{
        MakeModelObject({
            { { 0.0, 0.0, 0.0 }, Element::CARBON },
            { { 0.0, 0.0, 0.5 }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };

    const auto q_score{
        rhbm_gem::core::CalculateQScoreForAtom(
            atom,
            map,
            *model,
            height,
            offset,
            sigma,
            point_count)
    };

    constexpr double expected_q_score{ 0.29089961997690711 };
    EXPECT_NEAR(q_score, expected_q_score, 1.0e-6);
}

TEST(QScoreHelperTest, ReturnsZeroQScoreForConstantMap)
{
    constexpr double height{ 0.5 };
    constexpr double offset{ 1.0 };
    constexpr double sigma{ 0.5 };
    constexpr int point_count{ 4 };
    const auto map{
        MakeMapObject(
            { 5, 5, 5 },
            { 1.0, 1.0, 1.0 },
            { -2.0, -2.0, -2.0 },
            [](const std::array<double, 3> &)
            {
                return 3.5;
            })
    };
    const auto model{
        MakeModelObject({
            { { 0.0, 0.0, 0.0 }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };

    EXPECT_DOUBLE_EQ(
        rhbm_gem::core::CalculateQScoreForAtom(
            atom,
            map,
            *model,
            height,
            offset,
            sigma,
            point_count),
        0.0);
}

TEST(QScoreHelperTest, ReturnsZeroQScoreWhenEveryRadialShellIsEmpty)
{
    const auto map{
        MakeMapObject(
            { 5, 5, 5 },
            { 0.5, 0.5, 0.5 },
            { -1.0, -1.0, -1.0 },
            [](const std::array<double, 3> & position)
            {
                return position.at(2);
            })
    };
    const auto model{ MakeCrowdedModelObject() };
    const auto & atom{ *model->GetAtomList().front() };

    EXPECT_DOUBLE_EQ(
        rhbm_gem::core::CalculateQScoreForAtom(
            atom,
            map,
            *model,
            0.5,
            1.0,
            1.0,
            8),
        0.0);
}

TEST(QScoreHelperTest, RejectsSamplingPositionsOutsideMapBoundary)
{
    const auto map{
        MakeMapObject(
            { 3, 3, 3 },
            { 1.0, 1.0, 1.0 },
            { -1.0, -1.0, -1.0 },
            [](const std::array<double, 3> & position)
            {
                return position.at(2);
            })
    };
    const auto outside_center_model{
        MakeModelObject({
            { { 0.0, 0.0, 2.0 }, Element::CARBON }
        })
    };
    const auto radial_outside_model{
        MakeModelObject({
            { { 0.0, 0.0, 1.0 }, Element::CARBON }
        })
    };

    EXPECT_THROW(
        rhbm_gem::core::CalculateQScoreForAtom(
            *outside_center_model->GetAtomList().front(),
            map,
            *outside_center_model,
            0.5,
            0.5,
            0.5,
            2),
        std::out_of_range);
    EXPECT_THROW(
        rhbm_gem::core::CalculateQScoreForAtom(
            *radial_outside_model->GetAtomList().front(),
            map,
            *radial_outside_model,
            0.5,
            1.0,
            1.0,
            2),
        std::out_of_range);
}

TEST(QScoreHelperTest, RejectsInvalidQScoreArguments)
{
    const auto map{
        MakeMapObject(
            { 5, 5, 5 },
            { 0.5, 0.5, 0.5 },
            { -1.0, -1.0, -1.0 },
            [](const std::array<double, 3> & position)
            {
                return position.at(2);
            })
    };
    const auto model{
        MakeModelObject({
            { { 0.0, 0.0, 0.0 }, Element::CARBON }
        })
    };
    const auto other_model{
        MakeModelObject({
            { { 0.0, 0.0, 0.0 }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };
    const auto nan{ std::numeric_limits<double>::quiet_NaN() };
    const auto infinity{ std::numeric_limits<double>::infinity() };
    const std::array<double, 4> invalid_sigma_values{
        0.0,
        -1.0,
        nan,
        infinity
    };

    for (const auto invalid_sigma : invalid_sigma_values)
    {
        EXPECT_THROW(
            rhbm_gem::core::CalculateQScoreForAtom(
                atom, map, *model, 0.5, 1.0, invalid_sigma, 8),
            std::invalid_argument);
    }
    for (const auto invalid_num_points : { 1, 0, -1 })
    {
        EXPECT_THROW(
            rhbm_gem::core::CalculateQScoreForAtom(
                atom, map, *model, 0.5, 1.0, 0.1, invalid_num_points),
            std::invalid_argument);
    }
    EXPECT_THROW(
        rhbm_gem::core::CalculateQScoreForAtom(
            atom, map, *other_model, 0.5, 1.0, 0.1, 8),
        std::invalid_argument);
}

TEST(QScoreHelperTest, CalculatesSingleAtomAverageWithMapQParameters)
{
    const auto map{ MakeGaussianMapObject() };
    const auto model{
        MakeModelObject({
            { { 0.0, 0.0, 0.0 }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };
    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };
    const auto expected_q_score{
        rhbm_gem::core::CalculateQScoreForAtom(
            atom,
            map,
            *model,
            height,
            offset)
    };
    std::unordered_map<int, double> q_scores_by_serial_id;

    EXPECT_DOUBLE_EQ(
        rhbm_gem::core::CalculateAverageQScores(
            map,
            *model,
            height,
            offset,
            q_scores_by_serial_id),
        expected_q_score);
}

TEST(QScoreHelperTest, AveragesAllNonHydrogenAtomsRegardlessOfSelection)
{
    const auto map{ MakeGaussianMapObject() };
    const auto model{
        MakeModelObject({
            { { -0.5, 0.0, 0.0 }, Element::CARBON },
            { { 0.8, 0.0, 0.0 }, Element::OXYGEN },
            { { 100.0, 0.0, 0.0 }, Element::HYDROGEN }
        })
    };
    model->SelectAllAtoms();
    model->SetAtomSelected(2, false);
    ASSERT_EQ(model->GetSelectedAtomCount(), 2u);
    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    const auto first_q_score{
        rhbm_gem::core::CalculateQScoreForAtom(
            *model->GetAtomList().at(0),
            map,
            *model,
            height,
            offset)
    };
    const auto second_q_score{
        rhbm_gem::core::CalculateQScoreForAtom(
            *model->GetAtomList().at(1),
            map,
            *model,
            height,
            offset)
    };
    std::unordered_map<int, double> q_scores_by_serial_id;

    EXPECT_DOUBLE_EQ(
        rhbm_gem::core::CalculateAverageQScores(
            map,
            *model,
            height,
            offset,
            q_scores_by_serial_id),
        (first_q_score + second_q_score) / 2.0);
}

TEST(QScoreHelperTest, AverageOverloadUsesProvidedParametersAndReturnsEachNonHydrogenAtomScore)
{
    const auto map{ MakeGaussianMapObject() };
    const auto model{
        MakeModelObject({
            { { -0.5, 0.0, 0.0 }, Element::CARBON },
            { { 0.8, 0.0, 0.0 }, Element::OXYGEN },
            { { 0.0, 0.0, 0.0 }, Element::HYDROGEN }
        })
    };
    model->SelectAllAtoms(false);
    constexpr double height{ 1.75 };
    constexpr double offset{ -0.25 };
    const auto expected_first{
        rhbm_gem::core::CalculateQScoreForAtom(
            *model->GetAtomList().at(0), map, *model, height, offset)
    };
    const auto expected_second{
        rhbm_gem::core::CalculateQScoreForAtom(
            *model->GetAtomList().at(1), map, *model, height, offset)
    };
    std::unordered_map<int, double> q_scores_by_serial_id;

    const auto average{
        rhbm_gem::core::CalculateAverageQScores(
            map,
            *model,
            height,
            offset,
            q_scores_by_serial_id)
    };

    ASSERT_EQ(q_scores_by_serial_id.size(), 2u);
    EXPECT_DOUBLE_EQ(q_scores_by_serial_id.at(1), expected_first);
    EXPECT_DOUBLE_EQ(q_scores_by_serial_id.at(2), expected_second);
    EXPECT_EQ(q_scores_by_serial_id.count(3), 0u);
    EXPECT_DOUBLE_EQ(average, (expected_first + expected_second) / 2.0);
}

TEST(QScoreHelperTest, ReturnsZeroAverageWithoutNonHydrogenAtoms)
{
    const auto map{ MakeGaussianMapObject() };
    const rhbm_gem::ModelObject empty_model;
    const auto hydrogen_model{
        MakeModelObject({
            { { 100.0, 0.0, 0.0 }, Element::HYDROGEN }
        })
    };
    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    std::unordered_map<int, double> q_scores_by_serial_id{
        { 99, 1.0 }
    };
    EXPECT_DOUBLE_EQ(
        rhbm_gem::core::CalculateAverageQScores(
            map,
            empty_model,
            height,
            offset,
            q_scores_by_serial_id),
        0.0);
    EXPECT_TRUE(q_scores_by_serial_id.empty());

    q_scores_by_serial_id = { { 99, 1.0 } };
    EXPECT_DOUBLE_EQ(
        rhbm_gem::core::CalculateAverageQScores(
            map,
            *hydrogen_model,
            height,
            offset,
            q_scores_by_serial_id),
        0.0);
    EXPECT_TRUE(q_scores_by_serial_id.empty());
}

TEST(QScoreHelperTest, PropagatesOutOfRangeFromNonHydrogenAtoms)
{
    const auto map{
        MakeMapObject(
            { 5, 5, 5 },
            { 1.0, 1.0, 1.0 },
            { -2.0, -2.0, -2.0 },
            [](const std::array<double, 3> & position)
            {
                return position.at(2);
            })
    };
    const auto outside_center_model{
        MakeModelObject({
            { { 0.0, 0.0, 100.0 }, Element::CARBON }
        })
    };
    const auto radial_outside_model{
        MakeModelObject({
            { { 0.0, 0.0, 1.0 }, Element::CARBON }
        })
    };
    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    std::unordered_map<int, double> q_scores_by_serial_id{
        { 99, 1.0 }
    };
    EXPECT_THROW(
        rhbm_gem::core::CalculateAverageQScores(
            map,
            *outside_center_model,
            height,
            offset,
            q_scores_by_serial_id),
        std::out_of_range);
    EXPECT_EQ(
        q_scores_by_serial_id,
        (std::unordered_map<int, double>{ { 99, 1.0 } }));

    EXPECT_THROW(
        rhbm_gem::core::CalculateAverageQScores(
            map,
            *radial_outside_model,
            height,
            offset,
            q_scores_by_serial_id),
        std::out_of_range);
    EXPECT_EQ(
        q_scores_by_serial_id,
        (std::unordered_map<int, double>{ { 99, 1.0 } }));
}

TEST(QScoreHelperTest, RejectsTargetAtomFromDifferentModel)
{
    const auto target_model{
        MakeModelObject({
            { { 0.0, 0.0, 0.0 }, Element::CARBON }
        })
    };
    const auto other_model{
        MakeModelObject({
            { { 1.0, 0.0, 0.0 }, Element::CARBON }
        })
    };
    const auto & target_atom{ *target_model->GetAtomList().front() };

    EXPECT_THROW(
        rhbm_gem::core::GetRadialPointsForQScore(
            target_atom,
            *other_model,
            1.0,
            2),
        std::invalid_argument);
}

TEST(QScoreHelperTest, SharedInterpolationPreservesLegacyBoundaryDiagnostic)
{
    const auto map=MakeGaussianMapObject();
    const auto model=MakeModelObject({{{100.,0.,0.},Element::CARBON}});
    try
    {
        rhbm_gem::core::CalculateQScoreForAtom(*model->GetAtomList().front(),map,*model,1.,0.);
        FAIL()<<"Expected QScore boundary exception";
    }
    catch (const std::out_of_range & error)
    {
        EXPECT_STREQ(error.what(),"CalculateQScoreForAtom sampling position is outside the map boundary.");
    }
}

TEST(QScoreHelperTest, NeighborExclusionUsesStrictSquaredDistanceBoundary)
{
    // The north pole is at the origin, so its separation is exactly the neighbor's Z magnitude.
    constexpr double exclusion_radius{ 0.9 };
    const auto isolated{ MakeModelObject({ { { 0.0, 0.0, -1.0 }, Element::CARBON } }) };
    const auto reference{ rhbm_gem::core::GetRadialPointsForQScore(
        *isolated->GetAtomList().front(), *isolated, 1.0, 2) };
    ASSERT_EQ(reference.size(), 2u);
    ASSERT_EQ(reference.back().position, (std::array<double, 3>{}));

    for (const double separation : { std::nextafter(exclusion_radius, 0.0), exclusion_radius,
             std::nextafter(exclusion_radius, 1.0) })
    {
        const auto model{ MakeModelObject({
            { { 0.0, 0.0, -1.0 }, Element::CARBON },
            { { 0.0, 0.0, -separation }, Element::OXYGEN }
        }) };
        const auto shell{ rhbm_gem::core::GetRadialPointsForQScore(
            *model->GetAtomList().front(), *model, 1.0, 2) };
        if (separation < exclusion_radius)
        {
            EXPECT_GE(shell.size(), 2u);
            EXPECT_TRUE(std::none_of(shell.begin(), shell.end(), [&](const auto & point)
            {
                return point.position == reference.back().position;
            }));
        }
        else
        {
            ExpectSamplingPointListsEqual(shell, reference);
        }
    }
}

TEST(QScoreHelperTest, ZeroSquaredExclusionRadiusRetainsCoincidentPoint)
{
    constexpr double radius{ 1e-200 }; // Squaring the exclusion radius underflows to zero.
    const auto model{ MakeModelObject({
        { { 0.0, 0.0, 0.0 }, Element::CARBON },
        { { 0.0, 0.0, radius }, Element::OXYGEN }
    }) };
    const auto & atom{ *model->GetAtomList().front() };
    ASSERT_EQ(model->FindNeighborAtoms(atom, 2.0 * radius).size(), 1u);
    const auto shell{ rhbm_gem::core::GetRadialPointsForQScore(atom, *model, radius, 2) };
    ASSERT_EQ(shell.size(), 2u);
    EXPECT_EQ(shell.back().position, model->GetAtomList().back()->GetPositionRef());
}

TEST(QScoreHelperTest, MultipleOffGridNeighborsExcludeBothPoles)
{
    const std::array<double, 3> center{ 0.137, -0.219, 0.153 };
    const auto isolated{ MakeModelObject({ { center, Element::CARBON } }) };
    const auto model{ MakeModelObject({
        { center, Element::CARBON },
        { { center[0], center[1], center[2] - 0.5 }, Element::OXYGEN },
        { { center[0], center[1], center[2] + 0.5 }, Element::CARBON },
        { { center[0] + 1.95, center[1], center[2] }, Element::CARBON }
    }) };
    auto reference{ rhbm_gem::core::GetRadialPointsForQScore(
        *isolated->GetAtomList().front(), *isolated, 1.0, 4) };
    ASSERT_EQ(reference.size(), 4u);
    reference.pop_back();
    reference.erase(reference.begin());
    const auto shell{ rhbm_gem::core::GetRadialPointsForQScore(
        *model->GetAtomList().front(), *model, 1.0, 2) };
    ExpectSamplingPointListsEqual(shell, reference);
}

TEST(QScoreHelperTest, RejectsNonfiniteNeighborAndOverflowingSquaredExclusionRadius)
{
    const auto infinite_neighbor{ MakeModelObject({
        { { 0.0, 0.0, 0.0 }, Element::CARBON },
        { { 0.0, std::numeric_limits<double>::infinity(), 0.0 }, Element::OXYGEN }
    }) };
    const auto finite_neighbor{ MakeModelObject({
        { { 0.0, 0.0, 0.0 }, Element::CARBON },
        { { 0.0, 1.0, 0.0 }, Element::OXYGEN }
    }) };
    // Only the search-radius square overflows for the first case, exposing the nonfinite neighbor.
    EXPECT_THROW(rhbm_gem::core::GetRadialPointsForQScore(
        *infinite_neighbor->GetAtomList().front(), *infinite_neighbor, 8e153, 2), std::invalid_argument);
    EXPECT_THROW(rhbm_gem::core::GetRadialPointsForQScore(
        *finite_neighbor->GetAtomList().front(), *finite_neighbor, 1e155, 2), std::invalid_argument);
}

TEST(QScoreHelperTest, NonfiniteCenterIsRejectedByMapInterpolation)
{
    const auto model{ MakeModelObject({
        { { std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0 }, Element::CARBON }
    }) };
    const auto & atom{ *model->GetAtomList().front() };
    const auto shell{ rhbm_gem::core::GetRadialPointsForQScore(atom, *model, 1.0, 2) };
    ASSERT_EQ(shell.size(), 2u);
    for (const auto & point : shell)
    {
        EXPECT_TRUE(std::isnan(point.position[0]));
    }
    const auto map{ MakeGaussianMapObject() };
    EXPECT_THROW(rhbm_gem::core::CalculateQScoreForAtom(atom, map, *model, 1.0, 0.0), std::out_of_range);
}

// Frozen arithmetic from the pre-refactor Standard implementation. This guards
// bit-level compatibility independently of the new preparation/reference path.
namespace {
using rhbm_gem::AtomObject;
using rhbm_gem::MapObject;
using rhbm_gem::ModelObject;
using rhbm_gem::core::GetRadialPointsForQScore;
namespace numeric_validation=rhbm_gem::numeric_validation;
namespace detail=rhbm_gem::core::detail;
constexpr double kMaximumRadiusTolerance{0.01};
double InterpolateQScoreResponse(const MapObject & map, const std::array<double, 3> & position)
{
    try
    {
        return detail::InterpolateTrilinear(map, position);
    }
    catch (const std::out_of_range &)
    {
        throw std::out_of_range(
            "CalculateQScoreForAtom sampling position is outside the map boundary.");
    }
}

double LegacyMeanSubtractedCorrelation(
    const std::vector<double> & map_values,
    const std::vector<double> & reference_values)
{
    if (map_values.empty())
    {
        return 0.0;
    }

    double map_value_sum{ 0.0 };
    double reference_value_sum{ 0.0 };
    for (std::size_t i = 0; i < map_values.size(); ++i)
    {
        map_value_sum += map_values.at(i);
        reference_value_sum += reference_values.at(i);
    }
    const auto sample_count{ static_cast<double>(map_values.size()) };
    const auto map_value_mean{ map_value_sum / sample_count };
    const auto reference_value_mean{ reference_value_sum / sample_count };

    double numerator{ 0.0 };
    double map_value_norm_square{ 0.0 };
    double reference_value_norm_square{ 0.0 };
    for (std::size_t i = 0; i < map_values.size(); ++i)
    {
        const auto map_value_difference{
            map_values.at(i) - map_value_mean
        };
        const auto reference_value_difference{
            reference_values.at(i) - reference_value_mean
        };
        numerator += map_value_difference * reference_value_difference;
        map_value_norm_square += map_value_difference * map_value_difference;
        reference_value_norm_square += reference_value_difference * reference_value_difference;
    }

    const auto denominator{
        std::sqrt(map_value_norm_square) * std::sqrt(reference_value_norm_square)
    };
    if (denominator == 0.0 ||
        !std::isfinite(denominator) ||
        !std::isfinite(numerator))
    {
        return 0.0;
    }

    const auto q_score{ numerator / denominator };
    return std::isfinite(q_score) ? q_score : 0.0;
}

[[gnu::noinline]] double LegacyStandardScore(
    const AtomObject & atom,
    const MapObject & map,
    const ModelObject & model,
    double height,
    double offset,
    double sigma,
    int num_points)
{
    numeric_validation::RequireFinitePositive(sigma, "CalculateQScoreForAtom sigma");
    numeric_validation::RequireAtLeast(num_points, 2, "CalculateQScoreForAtom num_points");

    const auto center_map_value{
        InterpolateQScoreResponse(map, atom.GetPositionRef())
    };
    const auto center_reference_value{ height + offset };

    std::vector<double> map_values(static_cast<std::size_t>(num_points), center_map_value);
    std::vector<double> reference_values(static_cast<std::size_t>(num_points), center_reference_value);

    for (double radius = 0.1; radius < 2.0 + kMaximumRadiusTolerance; radius += 0.1)
    {
        const auto radial_points{
            GetRadialPointsForQScore(atom, model, radius, num_points)
        };
        if (radial_points.empty()) continue;
        const auto reference_value{
            height * std::exp(-0.5 * radius * radius / (sigma * sigma)) + offset
        };
        for (const auto & radial_point : radial_points)
        {
            map_values.emplace_back(InterpolateQScoreResponse(map, radial_point.position));
            reference_values.emplace_back(reference_value);
        }
    }

    return LegacyMeanSubtractedCorrelation(map_values, reference_values);
}


}

TEST(QScoreHelperTest, PreparedSamplesAndStandardArithmeticMatchLegacyExactly)
{
    const auto map=MakeGaussianMapObject();
    for (const double neighbor_z:{.5,1.2,4.0})
    {
        const auto model=MakeModelObject({{{0,0,0}},{{0,0,neighbor_z}}});
        const auto & atom=*model->GetAtomList().front();
        for (const int count:{2,8})
        {
            const auto samples=rhbm_gem::core::PrepareQScoreSamples(atom,map,*model,count);
            std::size_t k{};
            for (int center=0;center<count;++center,++k)
            {
                EXPECT_EQ(samples[k].point.position,atom.GetPositionRef());
                EXPECT_EQ(samples[k].point.distance,0);
                EXPECT_EQ(samples[k].response,1);
            }
            for (double radius=.1;radius<2.01;radius+=.1)
                for (const auto & point:rhbm_gem::core::GetRadialPointsForQScore(atom,*model,radius,count))
                {
                    ASSERT_LT(k,samples.size());
                    EXPECT_EQ(samples[k].point.position,point.position);
                    EXPECT_EQ(samples[k].point.distance,radius);
                    EXPECT_EQ(samples[k].response,rhbm_gem::core::detail::InterpolateTrilinear(map,point.position));
                    ++k;
                }
            EXPECT_EQ(k,samples.size());
            EXPECT_EQ(rhbm_gem::core::CalculateQScoreForAtom(atom,map,*model,1.7,-.31,.6,count),
                LegacyStandardScore(atom,map,*model,1.7,-.31,.6,count));
        }
    }
}

TEST(QScoreHelperTest, ExplicitCorrelationDistinguishesZeroFromUnavailable)
{
    using rhbm_gem::core::CalculateQScoreCorrelation;
    const auto zero=CalculateQScoreCorrelation({-1,0,1},{1,-2,1});
    ASSERT_TRUE(zero.value); EXPECT_EQ(*zero.value,0);
    EXPECT_EQ(CalculateQScoreCorrelation({.1,.1,.1},{0,1,2}).reason,"observation-zero-variance");
    EXPECT_EQ(CalculateQScoreCorrelation({0,1,2},{.1,.1,.1}).reason,"reference-zero-variance");
    EXPECT_EQ(CalculateQScoreCorrelation({},{}).reason,"sample-unavailable");
    EXPECT_EQ(CalculateQScoreCorrelation({0,std::numeric_limits<double>::infinity()},{0,1}).reason,"nonfinite-response");
}
