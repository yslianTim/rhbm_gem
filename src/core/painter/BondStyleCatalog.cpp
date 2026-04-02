#include "core/painter/BondStyleCatalog.hpp"

#include "data/detail/BondClassifier.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <vector>

namespace rhbm_gem {

namespace {

const std::vector<StyleSpec> kMainChainMemberStyleList
{
    { 633, 21, 25, 1, "N-C_{#alpha}" },
    { 881, 20, 24, 1, "C_{#alpha}-C" },
    { 418, 22, 26, 1, "C=O" },
    { 862, 23, 32, 1, "C-N" }
};

const StyleSpec kUnknownStyle{};

bool IsValidMainChainMemberID(size_t id)
{
    if (id >= BondClassifier::GetMainChainMemberCount())
    {
        Logger::Log(LogLevel::Error, "Invalid main chain member ID: " + std::to_string(id));
        return false;
    }
    return true;
}

} // namespace

StyleSpec BondStyleCatalog::GetMainChainMemberStyle(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return kUnknownStyle;
    return kMainChainMemberStyleList.at(id);
}

} // namespace rhbm_gem
