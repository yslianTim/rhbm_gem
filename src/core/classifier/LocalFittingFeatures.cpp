#include "core/classifier/LocalFittingFeatures.hpp"

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace rhbm_gem::core::classifier {

namespace {

using Wide = boost::multiprecision::cpp_bin_float_quad;
constexpr double kMissingFeature{ std::numeric_limits<double>::quiet_NaN() };

double OptionalFeatureValue(const std::optional<double> & value)
{
    return value.value_or(kMissingFeature);
}

double NarrowFeatureValue(const Wide & value)
{
    if (abs(value) > std::numeric_limits<double>::max()) return kMissingFeature;
    const auto narrowed{ static_cast<double>(value) };
    return value < 0 ? std::copysign(narrowed, -1.0) : narrowed;
}

} // namespace

JointModelResponseFeatures EvaluateJointModelResponseFeatures(const GaussianModel3D & point)
{
    const double amplitude{ point.GetAmplitude() };
    const double width{ point.GetWidth() };
    const double charge{ point.GetOffset() };
    if (!std::isfinite(amplitude) || amplitude < 0.0 ||
        !std::isfinite(width) || width <= 0.0 || !std::isfinite(charge))
    {
        return { kMissingFeature, kMissingFeature, kMissingFeature };
    }

    constexpr double radius{ GetPeelingRegionBounds(PeelingRegion::Signal).max_distance };
    using boost::math::constants::pi;
    using boost::math::constants::root_two;
    using boost::math::constants::root_two_div_pi;
    static const Wide four_pi{ 4 * pi<Wide>() };
    const Wide inverse_width{ Wide{1} / width };
    // Factor out the charge basis. The wide exponent and significand preserve
    // cancellation even when either separate center contribution exceeds double.
    const Wide center{ root_two_div_pi<Wide>() * inverse_width *
        (Wide{amplitude} * inverse_width * inverse_width / four_pi + charge) };

    Wide gaussian_mass, charge_volume;
    if (width >= radius)
    {
        const Wide x{ radius * inverse_width };
        const double y{ -0.5 * static_cast<double>(x * x) };
        double gaussian_sum{}, charge_sum{};
        // Integrate the exp and erf series term by term. For x <= 1 the
        // first omitted term is bounded by (1/2)^15 / (15! * 33) in each sum.
        for (int n = 14; n >= 0; --n)
        {
            const double step{ y / (n + 1) };
            gaussian_sum = 1.0 / (2 * n + 3) + step * gaussian_sum;
            charge_sum = 1.0 / ((2 * n + 1) * (2 * n + 3)) + step * charge_sum;
        }
        gaussian_mass = root_two_div_pi<Wide>() * x * x * x * gaussian_sum;
        charge_volume = four_pi * root_two_div_pi<Wide>() * (radius * radius) * x * charge_sum;
    }
    else if (width <= radius / 64.0)
    {
        // At x >= 64 the omitted exponential tails are below binary64
        // precision. The charge integral still has its algebraic width correction.
        gaussian_mass = 1;
        charge_volume = (four_pi / 2) * (radius * radius - Wide{width} * width);
    }
    else
    {
        const double x{ radius / width };
        const double erf{ std::erf(x / root_two<double>()) };
        const double exponential{ root_two_div_pi<double>() * std::exp(-0.5 * x * x) };
        gaussian_mass = erf - x * exponential;
        charge_volume = (four_pi / 2) *
            ((radius * radius - width * width) * erf + radius * width * exponential);
    }
    return {
        NarrowFeatureValue(center),
        NarrowFeatureValue(gaussian_mass * amplitude),
        NarrowFeatureValue(charge_volume * charge),
    };
}

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

        const auto response{ point ? EvaluateJointModelResponseFeatures(*point) :
            JointModelResponseFeatures{ missing, missing, missing } };
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
            response.center_height,
            response.gaussian_integral,
            response.charge_integral,
        };
        row.tail_peeling_ratio_unavailable_reason = std::move(tail_assessment.reason);
        rows.emplace_back(std::move(row));
    }
    return rows;
}

} // namespace rhbm_gem::core::classifier
