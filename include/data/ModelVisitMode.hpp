#pragma once

namespace rhbm_gem {

enum class ModelVisitMode
{
    AtomsThenSelf,
    BondsThenSelf,
    AtomsAndBondsThenSelf,
    SelfOnly
};

} // namespace rhbm_gem
