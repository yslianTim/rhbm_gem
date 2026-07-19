#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace {

struct AtomSpec
{
    std::array<float, 3> position;
    Element element{ Element::CARBON };
};

rhbm_gem::MapObject MakeMapObject(
    const std::array<int, 3> & grid_size,
    const std::array<float, 3> & grid_spacing,
    const std::array<float, 3> & origin,
    const std::vector<float> & values)
{
    auto map_values{ std::make_unique<float[]>(values.size()) };
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

rhbm_gem::MapObject MakeMapObject(const std::vector<float> & values)
{
    return MakeMapObject(
        {
            static_cast<int>(values.size()),
            1,
            1
        },
        { 1.0f, 1.0f, 1.0f },
        { 0.0f, 0.0f, 0.0f },
        values);
}

rhbm_gem::MapObject MakeMapObject(
    const std::array<int, 3> & grid_size,
    const std::array<float, 3> & grid_spacing,
    const std::array<float, 3> & origin,
    const std::function<float(const std::array<float, 3> &)> & value_at_position)
{
    std::vector<float> values;
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
                        origin.at(0) + static_cast<float>(x) * grid_spacing.at(0),
                        origin.at(1) + static_cast<float>(y) * grid_spacing.at(1),
                        origin.at(2) + static_cast<float>(z) * grid_spacing.at(2)
                    }));
            }
        }
    }
    return MakeMapObject(grid_size, grid_spacing, origin, values);
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
        { { 0.0f, 0.0f, 0.0f }, Element::CARBON }
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
                    std::sqrt(static_cast<float>(x * x + y * y + z * z))
                };
                atom_specs.emplace_back(AtomSpec{
                    {
                        static_cast<float>(x) / norm,
                        static_cast<float>(y) / norm,
                        static_cast<float>(z) / norm
                    },
                    Element::CARBON
                });
            }
        }
    }
    return MakeModelObject(atom_specs);
}

double ComputeDistanceSquare(
    const std::array<float, 3> & lhs,
    const std::array<float, 3> & rhs)
{
    double distance_square{ 0.0 };
    for (std::size_t axis = 0; axis < lhs.size(); ++axis)
    {
        const auto difference{
            static_cast<double>(lhs.at(axis)) -
                static_cast<double>(rhs.at(axis))
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
        EXPECT_FLOAT_EQ(lhs.at(i).distance, rhs.at(i).distance);
        EXPECT_EQ(lhs.at(i).position, rhs.at(i).position);
        EXPECT_EQ(lhs.at(i).is_selected, rhs.at(i).is_selected);
    }
}

} // namespace

TEST(QScoreHelperTest, ReturnsAmplitudeAndOffsetForTypicalMap)
{
    const auto map{ MakeMapObject({ 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f }) };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    constexpr float expected_standard_deviation{ 2.44948974f };
    constexpr float expected_offset{ 4.5f - expected_standard_deviation };
    constexpr float expected_height{ 8.0f - expected_offset };
    EXPECT_NEAR(height, expected_height, 1.0e-5f);
    EXPECT_NEAR(offset, expected_offset, 1.0e-5f);
}

TEST(QScoreHelperTest, CapsReferenceHighAtMeanPlusTenStandardDeviations)
{
    std::vector<float> values(121, 0.0f);
    values.back() = 1.0f;
    const auto map{ MakeMapObject(values) };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    constexpr float expected_height{ 1.0f / 121.0f + 10.0f / 11.0f };
    EXPECT_NEAR(height, expected_height, 1.0e-5f);
    EXPECT_FLOAT_EQ(offset, 0.0f);
    EXPECT_LT(height + offset, map.GetMapValueMax());
}

TEST(QScoreHelperTest, ClampsOffsetToObservedMinimum)
{
    const auto map{ MakeMapObject({ 0.0f, 0.0f, 0.0f, 10.0f }) };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    EXPECT_FLOAT_EQ(height, 10.0f);
    EXPECT_FLOAT_EQ(offset, 0.0f);
}

TEST(QScoreHelperTest, ReturnsZeroAmplitudeForConstantMap)
{
    const auto map{
        MakeMapObject({ 3.5f, 3.5f, 3.5f, 3.5f, 3.5f, 3.5f, 3.5f, 3.5f })
    };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    EXPECT_FLOAT_EQ(height, 0.0f);
    EXPECT_FLOAT_EQ(offset, 3.5f);
}

TEST(QScoreHelperTest, ReturnsDeterministicSpiralPointsForIsolatedAtom)
{
    const auto model{
        MakeModelObject({
            { { 1.0f, 2.0f, 3.0f }, Element::CARBON }
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
        EXPECT_FLOAT_EQ(point.distance, 1.0f);
        EXPECT_TRUE(point.is_selected);
        EXPECT_NEAR(
            ComputeDistanceSquare(point.position, atom.GetPositionRef()),
            1.0,
            1.0e-5);
    }
    EXPECT_NEAR(first_points.front().position.at(0), 1.0f, 1.0e-6f);
    EXPECT_NEAR(first_points.front().position.at(1), 2.0f, 1.0e-6f);
    EXPECT_NEAR(first_points.front().position.at(2), 2.0f, 1.0e-6f);
    EXPECT_NEAR(first_points.back().position.at(0), 1.0f, 1.0e-6f);
    EXPECT_NEAR(first_points.back().position.at(1), 2.0f, 1.0e-6f);
    EXPECT_NEAR(first_points.back().position.at(2), 4.0f, 1.0e-6f);
}

TEST(QScoreHelperTest, ReturnsAllAcceptedPointsFromFirstSuccessfulRetry)
{
    const auto model{
        MakeModelObject({
            { { 0.0f, 0.0f, 0.0f }, Element::CARBON },
            { { 0.0f, 0.0f, 0.5f }, Element::CARBON }
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
            { { 0.0f, 0.0f, 0.0f }, Element::CARBON },
            { { 0.0f, 0.0f, 0.5f }, Element::HYDROGEN }
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
    EXPECT_NEAR(points.front().position.at(2), -1.0f, 1.0e-6f);
    EXPECT_NEAR(points.back().position.at(2), 1.0f, 1.0e-6f);
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
            { { 0.0f, 0.0f, 0.0f }, Element::CARBON }
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
    constexpr double sigma{ 0.4 };
    const auto map{
        MakeMapObject(
            { 41, 41, 41 },
            { 0.1f, 0.1f, 0.1f },
            { -2.0f, -2.0f, -2.0f },
            [](const std::array<float, 3> & position)
            {
                const auto distance_square{
                    static_cast<double>(position.at(0)) *
                        static_cast<double>(position.at(0)) +
                    static_cast<double>(position.at(1)) *
                        static_cast<double>(position.at(1)) +
                    static_cast<double>(position.at(2)) *
                        static_cast<double>(position.at(2))
                };
                return static_cast<float>(
                    std::exp(-0.5 * distance_square / (sigma * sigma)));
            })
    };
    const auto model{
        MakeModelObject({
            { { 0.0f, 0.0f, 0.0f }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };

    const auto matching_score{
        rhbm_gem::core::CalculateQScoreForAtom(
            atom,
            map,
            *model,
            sigma,
            1.0,
            0.1,
            8)
    };
    const auto mismatched_score{
        rhbm_gem::core::CalculateQScoreForAtom(
            atom,
            map,
            *model,
            0.8,
            1.0,
            0.1,
            8)
    };

    EXPECT_GT(matching_score, 0.99);
    EXPECT_GT(matching_score, mismatched_score);
}

TEST(QScoreHelperTest, UsesCenterWeightMaximumShellAndAllAcceptedPoints)
{
    const auto map{
        MakeMapObject(
            { 5, 5, 5 },
            { 0.5f, 0.5f, 0.5f },
            { -1.0f, -1.0f, -1.0f },
            [](const std::array<float, 3> & position)
            {
                return position.at(2);
            })
    };
    const auto model{
        MakeModelObject({
            { { 0.0f, 0.0f, 0.0f }, Element::CARBON },
            { { 0.0f, 0.0f, 0.5f }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };

    const auto q_score{
        rhbm_gem::core::CalculateQScoreForAtom(
            atom,
            map,
            *model,
            0.5,
            1.0,
            1.0,
            2)
    };

    constexpr double expected_q_score{ 0.3611575593 };
    EXPECT_NEAR(q_score, expected_q_score, 1.0e-6);
}

TEST(QScoreHelperTest, ReturnsZeroQScoreForConstantMap)
{
    const auto map{
        MakeMapObject(
            { 5, 5, 5 },
            { 0.5f, 0.5f, 0.5f },
            { -1.0f, -1.0f, -1.0f },
            [](const std::array<float, 3> &)
            {
                return 3.5f;
            })
    };
    const auto model{
        MakeModelObject({
            { { 0.0f, 0.0f, 0.0f }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };

    EXPECT_DOUBLE_EQ(
        rhbm_gem::core::CalculateQScoreForAtom(
            atom,
            map,
            *model,
            0.5,
            1.0,
            0.5,
            4),
        0.0);
}

TEST(QScoreHelperTest, ReturnsZeroQScoreWhenEveryRadialShellIsEmpty)
{
    const auto map{
        MakeMapObject(
            { 5, 5, 5 },
            { 0.5f, 0.5f, 0.5f },
            { -1.0f, -1.0f, -1.0f },
            [](const std::array<float, 3> & position)
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
            { 1.0f, 1.0f, 1.0f },
            { -1.0f, -1.0f, -1.0f },
            [](const std::array<float, 3> & position)
            {
                return position.at(2);
            })
    };
    const auto outside_center_model{
        MakeModelObject({
            { { 0.0f, 0.0f, 2.0f }, Element::CARBON }
        })
    };
    const auto radial_outside_model{
        MakeModelObject({
            { { 0.0f, 0.0f, 1.0f }, Element::CARBON }
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
            { 0.5f, 0.5f, 0.5f },
            { -1.0f, -1.0f, -1.0f },
            [](const std::array<float, 3> & position)
            {
                return position.at(2);
            })
    };
    const auto model{
        MakeModelObject({
            { { 0.0f, 0.0f, 0.0f }, Element::CARBON }
        })
    };
    const auto other_model{
        MakeModelObject({
            { { 0.0f, 0.0f, 0.0f }, Element::CARBON }
        })
    };
    const auto & atom{ *model->GetAtomList().front() };
    const auto nan{ std::numeric_limits<double>::quiet_NaN() };
    const auto infinity{ std::numeric_limits<double>::infinity() };
    const std::array<double, 4> invalid_positive_values{
        0.0,
        -1.0,
        nan,
        infinity
    };

    for (const auto invalid_value : invalid_positive_values)
    {
        EXPECT_THROW(
            rhbm_gem::core::CalculateQScoreForAtom(
                atom, map, *model, invalid_value, 1.0, 0.1, 8),
            std::invalid_argument);
        EXPECT_THROW(
            rhbm_gem::core::CalculateQScoreForAtom(
                atom, map, *model, 0.5, invalid_value, 0.1, 8),
            std::invalid_argument);
        EXPECT_THROW(
            rhbm_gem::core::CalculateQScoreForAtom(
                atom, map, *model, 0.5, 1.0, invalid_value, 8),
            std::invalid_argument);
    }
    EXPECT_THROW(
        rhbm_gem::core::CalculateQScoreForAtom(
            atom, map, *model, 0.5, 1.0, 2.0, 8),
        std::invalid_argument);
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

TEST(QScoreHelperTest, RejectsTargetAtomFromDifferentModel)
{
    const auto target_model{
        MakeModelObject({
            { { 0.0f, 0.0f, 0.0f }, Element::CARBON }
        })
    };
    const auto other_model{
        MakeModelObject({
            { { 1.0f, 0.0f, 0.0f }, Element::CARBON }
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
