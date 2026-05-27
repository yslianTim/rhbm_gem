#pragma once

#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

#include "data/detail/AtomClassifier.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TAxis.h>
#endif

namespace rhbm_gem::painter_internal {

namespace style_detail {

struct StyleSpec
{
    short color{ 1 };
    short solid_marker{ 1 };
    short open_marker{ 1 };
    short line_style{ 1 };
    std::string label{ "UNK" };
};

inline const std::array<StyleSpec, 16> kMainChainElementStyleList{{
    { 633, 21, 25, 1, "C_{#alpha}" },
    { 881, 20, 24, 1, "C" },
    { 418, 22, 26, 1, "N" },
    { 862, 23, 32, 1, "O" },
    { 111, 33, 27, 1, "P" },
    { 862, 23, 32, 1, "O1" },
    { 862, 23, 32, 1, "O2" },
    { 862, 23, 32, 1, "O5'" },
    { 862, 23, 32, 1, "O4'" },
    { 862, 23, 32, 1, "O3'" },
    { 862, 23, 32, 1, "O2'" },
    { 881, 20, 24, 1, "C5'" },
    { 881, 20, 24, 1, "C4'" },
    { 881, 20, 24, 1, "C3'" },
    { 881, 20, 24, 1, "C2'" },
    { 881, 20, 24, 1, "C1'" }
}};

inline const StyleSpec kUnknownStyle{};

inline bool IsValidMainChainMemberID(size_t member_id)
{
    if (member_id >= data_internal::GetMainChainMemberCount())
    {
        Logger::Log(
            LogLevel::Error,
            "Invalid main chain member ID: " + std::to_string(member_id));
        return false;
    }
    return true;
}

inline StyleSpec GetMainChainElementStyle(size_t member_id)
{
    if (IsValidMainChainMemberID(member_id) == false) return kUnknownStyle;
    return kMainChainElementStyleList.at(member_id);
}

inline StyleSpec GetMainChainSpotStyle(Spot spot)
{
    switch (spot)
    {
    case Spot::CA:
        return StyleSpec{ 633, 21, 25, 1, "C_{#alpha}" };
    case Spot::C:
        return StyleSpec{ 881, 20, 24, 2, "C" };
    case Spot::N:
        return StyleSpec{ 418, 22, 26, 3, "N" };
    case Spot::O:
        return StyleSpec{ 862, 23, 32, 4, "O" };
    default:
        return kUnknownStyle;
    }
}

inline StyleSpec GetMainChainStructureStyle(Structure structure)
{
    switch (structure)
    {
    case Structure::FREE:
        return StyleSpec{ 633, 25, 25, 1, "N" };
    case Structure::HELX_P:
        return StyleSpec{ 418, 26, 26, 2, "#alpha" };
    case Structure::SHEET:
        return StyleSpec{ 862, 32, 32, 3, "#beta" };
    default:
        return kUnknownStyle;
    }
}

} // namespace style_detail

inline short GetMainChainElementColor(size_t member_id)
{
    return style_detail::GetMainChainElementStyle(member_id).color;
}

inline short GetMainChainSpotColor(Spot spot)
{
    return style_detail::GetMainChainSpotStyle(spot).color;
}

inline short GetMainChainStructureColor(Structure structure)
{
    return style_detail::GetMainChainStructureStyle(structure).color;
}

inline short GetMainChainElementSolidMarker(size_t member_id)
{
    return style_detail::GetMainChainElementStyle(member_id).solid_marker;
}

inline short GetMainChainSpotSolidMarker(Spot spot)
{
    return style_detail::GetMainChainSpotStyle(spot).solid_marker;
}

inline short GetMainChainStructureMarker(Structure structure)
{
    return style_detail::GetMainChainStructureStyle(structure).solid_marker;
}

inline short GetMainChainElementOpenMarker(size_t member_id)
{
    return style_detail::GetMainChainElementStyle(member_id).open_marker;
}

inline short GetMainChainSpotOpenMarker(Spot spot)
{
    return style_detail::GetMainChainSpotStyle(spot).open_marker;
}

inline short GetMainChainSpotLineStyle(Spot spot)
{
    return style_detail::GetMainChainSpotStyle(spot).line_style;
}

inline short GetMainChainStructureLineStyle(Structure structure)
{
    return style_detail::GetMainChainStructureStyle(structure).line_style;
}

inline std::string GetMainChainElementLabel(size_t member_id)
{
    return style_detail::GetMainChainElementStyle(member_id).label;
}

inline std::string GetMainChainSpotLabel(Spot spot)
{
    return style_detail::GetMainChainSpotStyle(spot).label;
}

inline std::string GetMainChainStructureLabel(Structure structure)
{
    return style_detail::GetMainChainStructureStyle(structure).label;
}

inline std::string BuildPainterOutputLabel(const ModelObject & model_object)
{
    auto is_simulation{ model_object.GetEmdID().find("Simulation") != std::string::npos };
    auto label{ model_object.GetPdbID() };
    if (is_simulation)
    {
        label += "_" + model_object.GetEmdID()
                 + "_bw" + string_helper::ToStringWithPrecision(model_object.GetResolution(), 2);
    }
    label += ".pdf";
    return label;
}

#ifdef HAVE_ROOT
inline void RemodelAxisLabels(
    ::TAxis * axis,
    const std::vector<std::string> & label_list,
    double angle,
    int align)
{
    auto label_size{ static_cast<int>(label_list.size()) };
    root_helper::SetAxisTickAttribute(axis, 0.0, label_size + 1);
    axis->SetLimits(-0.75, static_cast<float>(label_size) - 0.5f);
    for (size_t i = 0; i < label_list.size(); i++)
    {
        auto label_index{ static_cast<int>(i) + 1 };
        axis->ChangeLabel(label_index, angle, -1, align, -1, -1, label_list.at(i).data());
    }
}
#endif

} // namespace rhbm_gem::painter_internal
