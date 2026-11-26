#include "core/classifier/LocalFittingFeatures.hpp"

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace rhbm_gem::core::classifier {

namespace {

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

std::vector<LocalFittingFeatureRow> BuildLocalFittingFeatureRows(const ModelObject & model_object)
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
        if (!local_view.IsAvailable())
        {
            LocalFittingFeatureRow row; row.serial_id=atom->GetSerialID(); row.residue=ChemicalDataHelper::GetLabel(atom->GetResidue()); row.spot=atom->GetAtomID();
            row.features.fill(std::numeric_limits<double>::quiet_NaN()); rows.push_back(std::move(row)); continue;
        }
        const auto & joint = local_view.GetJointStageEstimate();
        if (joint.source.method == EstimateMethod::JointComponents && joint.source.role != FittingRole::OutputTarget) continue;
        const auto & point = joint.point;
        const double missing=std::numeric_limits<double>::quiet_NaN();

        const auto signal_peeling_ratio{
            local_view.AssessLocalFittingPeelingRatio(PeelingRegion::Signal).ratio
        };
        auto tail_assessment{ local_view.AssessLocalFittingPeelingRatio(PeelingRegion::Tail) };

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
        const double distance_to_closest_neighbor = closest_neighbors.empty() ? missing :
            array_helper::ComputeNorm(closest_neighbors.front()->GetPositionRef(), position);

        double neighbor_distance_sum{ 0.0 };
        for (const auto * neighbor : neighbors)
        {
            if (neighbor == atom) continue;
            const auto & neighbor_position{ neighbor->GetPositionRef() };
            const auto distance{
                array_helper::ComputeNorm(neighbor_position, position)
            };
            neighbor_distance_sum += distance;
        }

        LocalFittingFeatureRow row;
        row.serial_id = atom->GetSerialID();
        row.residue = ChemicalDataHelper::GetLabel(atom->GetResidue());
        row.spot = atom->GetAtomID();
        row.features = {
            neighbor_distance_sum,
            distance_to_closest_neighbor,
            OptionalFeatureValue(signal_peeling_ratio),
            OptionalFeatureValue(tail_assessment.ratio),
            point ? point->GetAmplitude() : missing,
            point ? point->GetWidth() : missing,
            point ? point->GetOffset() : missing,
        };
        row.tail_peeling_ratio_unavailable_reason = std::move(tail_assessment.reason);
        rows.emplace_back(std::move(row));
    }
    return rows;
}

} // namespace rhbm_gem::core::classifier
