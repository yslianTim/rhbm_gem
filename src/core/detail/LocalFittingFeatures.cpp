#include "core/detail/LocalFittingFeatures.hpp"

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include <algorithm>
#include <array>
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
    FittingStage stage,
    GaussianParameterGetter parameter_getter)
{
    const auto & current_model{
        AtomLocalPotentialView::For(atom).GetEstimateMDPDE(stage)
    };
    const auto current_value{ (current_model.*parameter_getter)() };
    int rank{ 1 };
    for (const auto * comparison_atom : comparison_atoms)
    {
        const auto & comparison_model{
            AtomLocalPotentialView::For(*comparison_atom).GetEstimateMDPDE(stage)
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

    auto kd_tree_root{ KDTreeAlgorithm<AtomObject>::BuildKDTree(atom_list) };
    std::sort(
        atom_list.begin(),
        atom_list.end(),
        [](const AtomObject * lhs, const AtomObject * rhs)
        {
            return lhs->GetSerialID() < rhs->GetSerialID();
        });

    constexpr std::array fitting_stages{
        FittingStage::First,
        FittingStage::Second,
        FittingStage::Third
    };

    std::vector<LocalFittingFeatureRow> rows;
    rows.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        const auto & first_model{ local_view.GetEstimateMDPDE(FittingStage::First) };
        const auto & second_model{ local_view.GetEstimateMDPDE(FittingStage::Second) };
        const auto & third_model{ local_view.GetEstimateMDPDE(FittingStage::Third) };

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

        std::array<int, fitting_stages.size()> amplitude_ranks{};
        std::array<int, fitting_stages.size()> width_ranks{};
        std::array<int, fitting_stages.size()> offset_ranks{};
        for (std::size_t stage_index = 0;
            stage_index < fitting_stages.size();
            ++stage_index)
        {
            const auto stage{ fitting_stages[stage_index] };
            amplitude_ranks[stage_index] = ComputeLocalParameterRank(
                *atom, comparison_atoms, stage, &GaussianModel3D::GetAmplitude);
            width_ranks[stage_index] = ComputeLocalParameterRank(
                *atom, comparison_atoms, stage, &GaussianModel3D::GetWidth);
            offset_ranks[stage_index] = ComputeLocalParameterRank(
                *atom, comparison_atoms, stage, &GaussianModel3D::GetOffset);
        }

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

        LocalFittingFeatureRow row;
        row.serial_id = atom->GetSerialID();
        row.residue = ChemicalDataHelper::GetLabel(atom->GetResidue());
        row.spot = atom->GetAtomID();
        row.features = {
            static_cast<double>(local_view.GetNeighborCountForPeeling()),
            static_cast<double>(atom->FindNeighborAtoms(2.0, false).size()),
            OptionalFeatureValue(signal_peeling_ratio),
            OptionalFeatureValue(tail_peeling_ratio),
            first_model.GetAmplitude(),
            second_model.GetAmplitude(),
            third_model.GetAmplitude(),
            first_model.GetWidth(),
            second_model.GetWidth(),
            third_model.GetWidth(),
            first_model.GetOffset(),
            second_model.GetOffset(),
            third_model.GetOffset(),
            static_cast<double>(amplitude_ranks[0]),
            static_cast<double>(amplitude_ranks[1]),
            static_cast<double>(amplitude_ranks[2]),
            static_cast<double>(width_ranks[0]),
            static_cast<double>(width_ranks[1]),
            static_cast<double>(width_ranks[2]),
            static_cast<double>(offset_ranks[0]),
            static_cast<double>(offset_ranks[1]),
            static_cast<double>(offset_ranks[2]),
        };
        rows.emplace_back(std::move(row));
    }
    return rows;
}

} // namespace rhbm_gem::core::detail
