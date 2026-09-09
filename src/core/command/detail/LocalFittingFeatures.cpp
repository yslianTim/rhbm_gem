#include "core/command/detail/LocalFittingFeatures.hpp"

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace rhbm_gem::core::detail {

namespace {

constexpr std::size_t kLocalRankNeighborCount{ 3 };
constexpr double kSignalPeelingDistanceMin{ 0.0 };
constexpr double kSignalPeelingDistanceMaxExclusive{ 1.0 };
constexpr double kTailPeelingDistanceMin{ 1.0 };
constexpr double kTailPeelingDistanceMax{ 2.0 };

using GaussianParameterGetter = double (GaussianModel3D::*)() const;

int ComputeLocalParameterRank(
    const AtomObject & atom,
    const std::vector<AtomObject *> & comparison_atoms,
    GaussianParameterGetter parameter_getter)
{
    const auto & current_model{
        AtomLocalPotentialView::For(atom).GetEstimateMDPDE(FittingStage::Second)
    };
    const auto current_value{ (current_model.*parameter_getter)() };
    int rank{ 1 };
    for (const auto * comparison_atom : comparison_atoms)
    {
        const auto & comparison_model{
            AtomLocalPotentialView::For(*comparison_atom).GetEstimateMDPDE(FittingStage::Second)
        };
        if ((comparison_model.*parameter_getter)() > current_value)
        {
            ++rank;
        }
    }
    return rank;
}

double OptionalFeatureValue(const std::optional<double> & value)
{
    return value.value_or(std::numeric_limits<double>::quiet_NaN());
}

} // namespace

std::string BuildLocalFittingCsvHeader()
{
    std::string header;
    const auto append_name = [&header](std::string_view name)
    {
        if (!header.empty()) header += ',';
        header += name;
    };
    for (const auto name : kLocalFittingIdentifierNames) append_name(name);
    for (const auto name : kLocalFittingFeatureNames) append_name(name);
    return header;
}

std::vector<LocalFittingFeatureRow> BuildLocalFittingFeatureRows(
    const ModelObject & model_object,
    bool peeling_applied)
{
    auto atom_list{ model_object.GetSelectedAtoms() };
    if (atom_list.empty()) return {};

    std::vector<AtomObject *> non_hydrogen_atoms;
    non_hydrogen_atoms.reserve(model_object.GetAtomList().size());
    for (const auto & atom : model_object.GetAtomList())
    {
        if (atom->GetElement() != Element::HYDROGEN)
        {
            non_hydrogen_atoms.push_back(atom.get());
        }
    }
    auto non_hydrogen_kd_tree_root{
        KDTreeAlgorithm<AtomObject>::BuildKDTree(non_hydrogen_atoms)
    };

    auto kd_tree_root{ KDTreeAlgorithm<AtomObject>::BuildKDTree(atom_list) };
    std::sort(
        atom_list.begin(),
        atom_list.end(),
        [](const AtomObject * lhs, const AtomObject * rhs)
        {
            return lhs->GetSerialID() < rhs->GetSerialID();
        });

    std::vector<LocalFittingFeatureRow> rows;
    rows.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        const auto & second_model{ local_view.GetEstimateMDPDE(FittingStage::Second) };

        auto comparison_atoms{ KDTreeAlgorithm<AtomObject>::KNearestNeighbors(
            kd_tree_root.get(),
            atom,
            std::min(kLocalRankNeighborCount + 1, atom_list.size()))
        };
        comparison_atoms.erase(
            std::remove(comparison_atoms.begin(), comparison_atoms.end(), atom),
            comparison_atoms.end());
        if (comparison_atoms.size() > kLocalRankNeighborCount)
        {
            comparison_atoms.resize(kLocalRankNeighborCount);
        }
        comparison_atoms.emplace_back(atom);

        const auto amplitude_rank{ ComputeLocalParameterRank(
            *atom, comparison_atoms, &GaussianModel3D::GetAmplitude) };
        const auto width_rank{ ComputeLocalParameterRank(
            *atom, comparison_atoms, &GaussianModel3D::GetWidth) };
        const auto offset_rank{ ComputeLocalParameterRank(
            *atom, comparison_atoms, &GaussianModel3D::GetOffset) };

        const auto signal_peeling_ratio{ local_view.GetLocalFittingPeelingRatio(
            peeling_applied,
            kSignalPeelingDistanceMin,
            std::nextafter(
                kSignalPeelingDistanceMaxExclusive,
                kSignalPeelingDistanceMin)) };
        const auto tail_peeling_ratio{ local_view.GetLocalFittingPeelingRatio(
            peeling_applied,
            kTailPeelingDistanceMin,
            kTailPeelingDistanceMax) };

        const auto neighbors{ KDTreeAlgorithm<AtomObject>::RangeSearch(
            non_hydrogen_kd_tree_root.get(), atom, 2.0)
        };
        const auto & position{ atom->GetPositionRef() };
        auto closest_neighbors{ KDTreeAlgorithm<AtomObject>::KNearestNeighbors(
            non_hydrogen_kd_tree_root.get(),
            atom,
            std::min(std::size_t{ 2 }, non_hydrogen_atoms.size()))
        };
        closest_neighbors.erase(
            std::remove(closest_neighbors.begin(), closest_neighbors.end(), atom),
            closest_neighbors.end());
        const auto & closest_position{ closest_neighbors.at(0)->GetPositionRef() };
        const auto distance_to_closest_neighbor{
            array_helper::ComputeNorm(closest_position, position)
        };

        std::size_t neighbor_count_in_2A{ 0 };
        std::size_t neighbor_count_in_1_5A{ 0 };
        double neighbor_distance_sum{ 0.0 };
        double neighbor_distance_sum_in_1_5A{ 0.0 };
        for (const auto * neighbor : neighbors)
        {
            if (neighbor == atom) continue;
            ++neighbor_count_in_2A;
            const auto & neighbor_position{ neighbor->GetPositionRef() };
            const auto distance{
                array_helper::ComputeNorm(neighbor_position, position)
            };
            neighbor_distance_sum += distance;
            if (distance <= 1.5)
            {
                ++neighbor_count_in_1_5A;
                neighbor_distance_sum_in_1_5A += distance;
            }
        }

        LocalFittingFeatureRow row;
        row.serial_id = atom->GetSerialID();
        row.residue = ChemicalDataHelper::GetLabel(atom->GetResidue());
        row.spot = atom->GetAtomID();
        row.features = {
            static_cast<double>(local_view.GetNeighborCountForPeeling()),
            static_cast<double>(neighbor_count_in_2A),
            static_cast<double>(neighbor_count_in_1_5A),
            neighbor_distance_sum,
            neighbor_distance_sum_in_1_5A,
            OptionalFeatureValue(signal_peeling_ratio),
            OptionalFeatureValue(tail_peeling_ratio),
            second_model.GetAmplitude(),
            second_model.GetWidth(),
            second_model.GetOffset(),
            static_cast<double>(amplitude_rank),
            static_cast<double>(width_rank),
            static_cast<double>(offset_rank),
            distance_to_closest_neighbor,
        };
        rows.emplace_back(std::move(row));
    }
    return rows;
}

} // namespace rhbm_gem::core::detail
