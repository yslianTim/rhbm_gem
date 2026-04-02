#pragma once

#include <cstddef>

#include "core/painter/StyleSpec.hpp"

namespace rhbm_gem {

class BondStyleCatalog
{
public:
    static StyleSpec GetMainChainMemberStyle(size_t id);

    static short GetMainChainMemberColor(size_t id)
    {
        return GetMainChainMemberStyle(id).color;
    }
    static short GetMainChainMemberSolidMarker(size_t id)
    {
        return GetMainChainMemberStyle(id).solid_marker;
    }
    static short GetMainChainMemberOpenMarker(size_t id)
    {
        return GetMainChainMemberStyle(id).open_marker;
    }
    static std::string GetMainChainMemberLabel(size_t id)
    {
        return GetMainChainMemberStyle(id).label;
    }
};

} // namespace rhbm_gem
