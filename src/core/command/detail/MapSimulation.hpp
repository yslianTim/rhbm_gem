#pragma once

#include <array>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <rhbm_gem/core/CommandTypes.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>

namespace rhbm_gem {
class MapObject;
class ModelObject;
}

namespace rhbm_gem::core::simulation {

struct SimulationAtom
{
    int serial_id{};
    int sequence_id{};
    std::string chain_id;
    std::string component_id;
    std::string atom_id;
    std::string alternate_indicator;
    Element element{ Element::UNK };
    Residue residue{ Residue::UNK };
    Spot spot{ Spot::UNK };
    Structure structure{ Structure::UNK };
    std::array<double, 3> position{};
    // No lookup is performed in neutral mode.
    std::optional<ChargeLookupResult> charge_lookup{};
    double charge_used{};
};

struct SimulationAtomPreparationResult
{
    std::vector<SimulationAtom> atom_list;
    std::array<double, 3> range_min{
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max() };
    std::array<double, 3> range_max{
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest() };
};

std::string_view ChargeStatusText(const SimulationAtom & atom);
SimulationAtomPreparationResult PrepareSimulationAtomList(
    const ModelObject & model, const MapSimulationRequest & request);
std::unique_ptr<MapObject> CreateMapObject(
    const MapSimulationRequest & request, const SimulationAtomPreparationResult & atoms);
// Returns the actual team size. Each z plane owns its voxels and sums atoms in preparation order.
int PopulateMapValueArray(MapObject & map, const SimulationAtomPreparationResult & atoms,
    const MapSimulationRequest & request, double blurring_width);

} // namespace rhbm_gem::core::simulation
