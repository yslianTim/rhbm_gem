#include "data/object/BondStyleCatalog.hpp"

#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <vector>

namespace rhbm_gem {

namespace {

const std::vector<short> kMainChainMemberColorList
{
    633, 881, 418, 862
};

const std::vector<short> kMainChainMemberSolidMarkerList
{
    21, 20, 22, 23
};

const std::vector<short> kMainChainMemberOpenMarkerList
{
    25, 24, 26, 32
};

const std::vector<std::string> kMainChainMemberTitleList
{
    "N-CA Bond",
    "CA-C Bond",
    "C=O Bond",
    "C-N Bond"
};

const std::vector<std::string> kMainChainMemberLabelList
{
    "N-C_{#alpha}",
    "C_{#alpha}-C",
    "C=O",
    "C-N"
};

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

short BondStyleCatalog::GetMainChainMemberColor(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return kMainChainMemberColorList.at(id);
}

short BondStyleCatalog::GetMainChainMemberSolidMarker(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return kMainChainMemberSolidMarkerList.at(id);
}

short BondStyleCatalog::GetMainChainMemberOpenMarker(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return kMainChainMemberOpenMarkerList.at(id);
}

const std::string & BondStyleCatalog::GetMainChainMemberLabel(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return kMainChainMemberLabelList.at(0);
    return kMainChainMemberLabelList.at(id);
}

const std::string & BondStyleCatalog::GetMainChainMemberTitle(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return kMainChainMemberTitleList.at(0);
    return kMainChainMemberTitleList.at(id);
}

} // namespace rhbm_gem
