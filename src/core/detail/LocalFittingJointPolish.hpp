#pragma once

#include "core/detail/LocalFittingTransformedChange.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

constexpr std::size_t kLocalFittingJointPolishShapeParameterSize{ 2 };

struct LocalFittingJointPolishParameterization
{
    std::vector<std::size_t> group_position_by_atom{};
    std::vector<std::vector<std::size_t>> atom_position_list_by_group{};
    Eigen::VectorXd seed_parameter{};

    std::size_t AtomCount() const
    {
        return group_position_by_atom.size();
    }

    std::size_t GroupCount() const
    {
        return atom_position_list_by_group.size();
    }

    Eigen::Index ParameterCount() const
    {
        return seed_parameter.size();
    }

    Eigen::Index ShapeColumn(
        std::size_t atom_position,
        std::size_t shape_parameter_index) const
    {
        return static_cast<Eigen::Index>(
            atom_position * kLocalFittingJointPolishShapeParameterSize +
            shape_parameter_index);
    }

    Eigen::Index OffsetColumn(std::size_t atom_position) const
    {
        return static_cast<Eigen::Index>(
            AtomCount() * kLocalFittingJointPolishShapeParameterSize +
            group_position_by_atom.at(atom_position));
    }

    std::optional<std::vector<GaussianModel3D>> DecodeModels(
        const Eigen::VectorXd & direction,
        double damping) const
    {
        if (direction.size() != ParameterCount() ||
            !std::isfinite(damping) || damping < 0.0 || damping > 1.0)
        {
            return std::nullopt;
        }
        const Eigen::VectorXd parameter{
            seed_parameter + damping * direction
        };
        if (!parameter.allFinite()) return std::nullopt;

        std::vector<GaussianModel3D> model_list;
        model_list.reserve(AtomCount());
        for (std::size_t atom_position = 0;
            atom_position < AtomCount();
            atom_position++)
        {
            const Eigen::Vector3d shape_coordinates{
                parameter(ShapeColumn(atom_position, 0)),
                parameter(ShapeColumn(atom_position, 1)),
                0.0
            };
            const auto shape_model{
                DecodeLocalFittingTransformedCoordinates(shape_coordinates)
            };
            if (!shape_model.has_value()) return std::nullopt;
            const auto model{
                shape_model->WithOffset(parameter(OffsetColumn(atom_position)))
            };
            if (!EncodeLocalFittingTransformedCoordinates(model).has_value())
            {
                return std::nullopt;
            }
            model_list.emplace_back(model);
        }
        return model_list;
    }
};

inline std::optional<LocalFittingJointPolishParameterization>
BuildLocalFittingJointPolishParameterization(
    const std::vector<GroupKey> & group_key_by_atom_position,
    const std::vector<GaussianModel3D> & base_model_list)
{
    if (group_key_by_atom_position.empty() ||
        group_key_by_atom_position.size() != base_model_list.size())
    {
        return std::nullopt;
    }

    std::map<GroupKey, std::size_t> group_position_by_key;
    for (const auto group_key : group_key_by_atom_position)
    {
        group_position_by_key.emplace(group_key, 0);
    }
    std::size_t group_position{ 0 };
    for (auto & [group_key, position] : group_position_by_key)
    {
        static_cast<void>(group_key);
        position = group_position++;
    }

    LocalFittingJointPolishParameterization parameterization;
    parameterization.group_position_by_atom.resize(base_model_list.size());
    parameterization.atom_position_list_by_group.resize(
        group_position_by_key.size());
    parameterization.seed_parameter = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(
            base_model_list.size() *
                kLocalFittingJointPolishShapeParameterSize +
            group_position_by_key.size()));
    std::vector<std::vector<double>> offset_list_by_group(
        group_position_by_key.size());

    for (std::size_t atom_position = 0;
        atom_position < base_model_list.size();
        atom_position++)
    {
        const auto transformed{
            EncodeLocalFittingTransformedCoordinates(
                base_model_list.at(atom_position))
        };
        if (!transformed.has_value()) return std::nullopt;
        const auto atom_group_position{
            group_position_by_key.at(group_key_by_atom_position.at(atom_position))
        };
        parameterization.group_position_by_atom.at(atom_position) =
            atom_group_position;
        parameterization.atom_position_list_by_group.at(atom_group_position)
            .emplace_back(atom_position);
        parameterization.seed_parameter(
            parameterization.ShapeColumn(atom_position, 0)) =
            (*transformed)(static_cast<Eigen::Index>(
                kLogPeakHeightChangeIndex));
        parameterization.seed_parameter(
            parameterization.ShapeColumn(atom_position, 1)) =
            (*transformed)(static_cast<Eigen::Index>(
                kLogWidthChangeIndex));
        offset_list_by_group.at(atom_group_position).emplace_back(
            base_model_list.at(atom_position).GetOffset());
    }

    for (std::size_t current_group_position = 0;
        current_group_position < offset_list_by_group.size();
        current_group_position++)
    {
        auto & offset_list{ offset_list_by_group.at(current_group_position) };
        if (offset_list.empty()) return std::nullopt;
        std::sort(offset_list.begin(), offset_list.end());
        const auto middle{ offset_list.size() / 2 };
        const auto median{
            offset_list.size() % 2 == 0 ?
                0.5 * offset_list.at(middle - 1) +
                    0.5 * offset_list.at(middle) :
                offset_list.at(middle)
        };
        if (!std::isfinite(median)) return std::nullopt;
        const auto representative_atom_position{
            parameterization.atom_position_list_by_group
                .at(current_group_position).front()
        };
        parameterization.seed_parameter(
            parameterization.OffsetColumn(representative_atom_position)) = median;
    }
    if (!parameterization.seed_parameter.allFinite()) return std::nullopt;
    return parameterization;
}

struct LocalFittingSharedOffsetResponse
{
    double response{ 0.0 };
    Eigen::Vector2d shape_jacobian{ Eigen::Vector2d::Zero() };
    double offset_jacobian{ 0.0 };
};

inline std::optional<LocalFittingSharedOffsetResponse>
EvaluateLocalFittingSharedOffsetResponse(
    const GaussianModel3D & model,
    double distance)
{
    if (!EncodeLocalFittingTransformedCoordinates(model).has_value() ||
        !std::isfinite(distance) || distance < 0.0)
    {
        return std::nullopt;
    }

    const auto width{ model.GetWidth() };
    const auto signal{ model.SignalAtDistance(distance) };
    const auto offset_basis{ model.OffsetBasisAtDistance(distance) };
    const auto response{ signal + model.GetOffset() * offset_basis };
    if (!std::isfinite(signal) || !std::isfinite(offset_basis) ||
        !std::isfinite(response))
    {
        return std::nullopt;
    }

    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    const auto normalized_distance{ distance / width };
    auto log_width_derivative{
        signal * normalized_distance * normalized_distance
    };
    if (distance < 1.0e-5)
    {
        log_width_derivative -=
            model.GetOffset() * center_offset_basis_scale / width;
    }
    else
    {
        const auto exponent{
            -0.5 * normalized_distance * normalized_distance
        };
        log_width_derivative -= model.GetOffset() *
            center_offset_basis_scale / width * std::exp(exponent);
    }

    Eigen::Vector2d shape_jacobian{
        signal,
        log_width_derivative
    };
    if (!shape_jacobian.allFinite()) return std::nullopt;

    return LocalFittingSharedOffsetResponse{
        response,
        shape_jacobian,
        offset_basis
    };
}

struct LocalFittingTransformedResponse
{
    double response{ 0.0 };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
};

struct LocalFittingTransformedModelInvariants
{
    GaussianModel3D model{};
    Eigen::Vector3d transformed{ Eigen::Vector3d::Zero() };
    double peak_height{ 0.0 };
};

inline std::optional<LocalFittingTransformedModelInvariants>
BuildLocalFittingTransformedModelInvariants(
    const GaussianModel3D & model)
{
    const auto transformed{ EncodeLocalFittingTransformedCoordinates(model) };
    if (!transformed.has_value()) return std::nullopt;

    const auto peak_height{
        std::exp((*transformed)(static_cast<Eigen::Index>(
            kLogPeakHeightChangeIndex)))
    };
    if (!std::isfinite(peak_height)) return std::nullopt;

    return LocalFittingTransformedModelInvariants{
        model,
        *transformed,
        peak_height
    };
}

inline std::optional<LocalFittingSharedOffsetResponse>
EvaluateLocalFittingSharedOffsetResponse(
    const LocalFittingTransformedModelInvariants & invariants,
    double distance)
{
    if (!std::isfinite(distance) || distance < 0.0) return std::nullopt;

    const auto evaluation{ invariants.model.EvaluateAtDistance(distance) };
    if (!std::isfinite(evaluation.signal) ||
        !std::isfinite(evaluation.offset_basis) ||
        !std::isfinite(evaluation.response))
    {
        return std::nullopt;
    }

    const auto width{ invariants.model.GetWidth() };
    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    const auto normalized_distance{ distance / width };
    auto log_width_derivative{
        evaluation.signal * normalized_distance * normalized_distance
    };
    if (distance < 1.0e-5)
    {
        log_width_derivative -=
            invariants.model.GetOffset() * center_offset_basis_scale / width;
    }
    else
    {
        const auto exponent{
            -0.5 * normalized_distance * normalized_distance
        };
        log_width_derivative -= invariants.model.GetOffset() *
            center_offset_basis_scale / width * std::exp(exponent);
    }

    const Eigen::Vector2d shape_jacobian{
        evaluation.signal,
        log_width_derivative
    };
    if (!shape_jacobian.allFinite()) return std::nullopt;

    return LocalFittingSharedOffsetResponse{
        evaluation.response,
        shape_jacobian,
        evaluation.offset_basis
    };
}

inline std::optional<LocalFittingTransformedResponse>
EvaluateLocalFittingTransformedResponse(
    const LocalFittingTransformedModelInvariants & invariants,
    double distance)
{
    const auto shared_offset_evaluation{
        EvaluateLocalFittingSharedOffsetResponse(invariants, distance)
    };
    if (!shared_offset_evaluation.has_value())
    {
        return std::nullopt;
    }

    const auto & model{ invariants.model };
    const auto width{ model.GetWidth() };
    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
    jacobian(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) =
        shared_offset_evaluation->shape_jacobian(0) +
        model.GetOffset() * shared_offset_evaluation->offset_jacobian;
    jacobian(static_cast<Eigen::Index>(kLogWidthChangeIndex)) =
        shared_offset_evaluation->shape_jacobian(1) +
        model.GetOffset() * shared_offset_evaluation->offset_jacobian;
    jacobian(static_cast<Eigen::Index>(kOffsetToPeakRatioChangeIndex)) =
        invariants.peak_height * width *
        shared_offset_evaluation->offset_jacobian /
        center_offset_basis_scale;
    if (!jacobian.allFinite())
    {
        return std::nullopt;
    }

    return LocalFittingTransformedResponse{
        shared_offset_evaluation->response,
        jacobian
    };
}

inline std::optional<LocalFittingTransformedResponse>
EvaluateLocalFittingTransformedResponse(
    const GaussianModel3D & model,
    double distance)
{
    const auto invariants{ BuildLocalFittingTransformedModelInvariants(model) };
    if (!invariants.has_value()) return std::nullopt;
    return EvaluateLocalFittingTransformedResponse(*invariants, distance);
}

} // namespace rhbm_gem::core::detail
