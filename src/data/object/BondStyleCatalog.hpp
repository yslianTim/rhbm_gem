#pragma once

#include <cstddef>
#include <string>

namespace rhbm_gem {

class BondStyleCatalog
{
public:
    static short GetMainChainMemberColor(size_t id);
    static short GetMainChainMemberSolidMarker(size_t id);
    static short GetMainChainMemberOpenMarker(size_t id);
    static const std::string & GetMainChainMemberLabel(size_t id);
    static const std::string & GetMainChainMemberTitle(size_t id);
};

} // namespace rhbm_gem
