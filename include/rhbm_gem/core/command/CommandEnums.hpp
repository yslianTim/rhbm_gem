#pragma once

namespace rhbm_gem {

enum class PainterType : int
{
    GAUS       = 0,
    MODEL      = 1,
    COMPARISON = 2,
    DEMO       = 3,
    ATOM       = 4
};

enum class PrinterType : int
{
    ATOM_POSITION  = 0,
    MAP_VALUE      = 1,
    GAUS_ESTIMATES = 2,
    ATOM_OUTLIER   = 3
};

enum class PotentialModel : int
{
    SINGLE_GAUS      = 0,
    FIVE_GAUS_CHARGE = 1,
    SINGLE_GAUS_USER = 2
};

enum class PartialCharge : int
{
    NEUTRAL = 0,
    PARTIAL = 1,
    AMBER   = 2
};

enum class TesterType : int
{
    BENCHMARK          = 0,
    DATA_OUTLIER       = 1,
    MEMBER_OUTLIER     = 2,
    MODEL_ALPHA_DATA   = 3,
    MODEL_ALPHA_MEMBER = 4,
    NEIGHBOR_DISTANCE  = 5
};

} // namespace rhbm_gem
