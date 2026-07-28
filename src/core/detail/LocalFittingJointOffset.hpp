#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

struct LocalFittingJointOffsetParameterization
{
    std::vector<std::size_t> group_position_by_atom{};
    std::vector<std::vector<std::size_t>> atom_position_list_by_group{};
    Eigen::VectorXd seed_offset{};

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
        return seed_offset.size();
    }

    Eigen::Index OffsetColumn(std::size_t atom_position) const
    {
        return static_cast<Eigen::Index>(
            group_position_by_atom.at(atom_position));
    }

    std::optional<Eigen::VectorXd> AggregateBasis(
        const std::vector<std::pair<std::size_t, double>> & atom_basis_entries) const
    {
        Eigen::VectorXd group_basis{ Eigen::VectorXd::Zero(ParameterCount()) };
        for (const auto & [atom_position, basis] : atom_basis_entries)
        {
            if (atom_position >= AtomCount() || !std::isfinite(basis))
            {
                return std::nullopt;
            }
            group_basis(OffsetColumn(atom_position)) += basis;
        }
        return group_basis.allFinite() ?
            std::optional<Eigen::VectorXd>{ std::move(group_basis) } :
            std::nullopt;
    }

    std::optional<Eigen::VectorXd> ExpandOffsets(
        const Eigen::VectorXd & group_offset) const
    {
        if (group_offset.size() != ParameterCount() || !group_offset.allFinite())
        {
            return std::nullopt;
        }

        Eigen::VectorXd atom_offset{
            Eigen::VectorXd::Zero(static_cast<Eigen::Index>(AtomCount()))
        };
        for (std::size_t atom_position = 0;
            atom_position < AtomCount();
            atom_position++)
        {
            atom_offset(static_cast<Eigen::Index>(atom_position)) =
                group_offset(OffsetColumn(atom_position));
        }
        return atom_offset;
    }
};

inline std::optional<LocalFittingJointOffsetParameterization>
BuildLocalFittingJointOffsetParameterization(
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

    LocalFittingJointOffsetParameterization parameterization;
    parameterization.group_position_by_atom.resize(base_model_list.size());
    parameterization.atom_position_list_by_group.resize(
        group_position_by_key.size());
    parameterization.seed_offset = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(group_position_by_key.size()));
    std::vector<std::vector<double>> offset_list_by_group(
        group_position_by_key.size());

    for (std::size_t atom_position = 0;
        atom_position < base_model_list.size();
        atom_position++)
    {
        const auto offset{ base_model_list.at(atom_position).GetOffset() };
        if (!std::isfinite(offset)) return std::nullopt;

        const auto atom_group_position{
            group_position_by_key.at(group_key_by_atom_position.at(atom_position))
        };
        parameterization.group_position_by_atom.at(atom_position) =
            atom_group_position;
        parameterization.atom_position_list_by_group.at(atom_group_position)
            .emplace_back(atom_position);
        offset_list_by_group.at(atom_group_position).emplace_back(offset);
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
        parameterization.seed_offset(
            static_cast<Eigen::Index>(current_group_position)) = median;
    }
    return parameterization.seed_offset.allFinite() ?
        std::optional<LocalFittingJointOffsetParameterization>{
            std::move(parameterization) } :
        std::nullopt;
}

} // namespace rhbm_gem::core::detail
