#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>

#ifndef UINT16_MAX
#define UINT16_MAX 65535
#endif

#ifndef UINT8_MAX
#define UINT8_MAX 255
#endif

enum class Residue : uint16_t
{
    ALA =  0, ARG =  1, ASN =  2, ASP =  3, CYS =  4,
    GLN =  5, GLU =  6, GLY =  7, HIS =  8, ILE =  9,
    LEU = 10, LYS = 11, MET = 12, PHE = 13, PRO = 14,
    SER = 15, THR = 16, TRP = 17, TYR = 18, VAL = 19,
    CSX = 20,
    HOH = 99,
    NA = 911, MG = 912, CL = 917, CA = 920, FE = 926, FE2 = 826, ZN = 930,
    UNK = UINT16_MAX
};

enum class Element : uint16_t
{
    HYDROGEN   =  1, HELIUM     =  2, LITHIUM    =  3, BERYLLIUM    =  4, BORON      =  5,
    CARBON     =  6, NITROGEN   =  7, OXYGEN     =  8, FLUORINE     =  9, NEON       = 10,
    SODIUM     = 11, MAGNESIUM  = 12, ALUMINUM   = 13, SILICON      = 14, PHOSPHORUS = 15,
    SULFUR     = 16, CHLORINE   = 17, ARGON      = 18, POTASSIUM    = 19, CALCIUM    = 20,
    SCANDIUM   = 21, TITANIUM   = 22, VANADIUM   = 23, CHROMIUM     = 24, MANGANESE  = 25,
    IRON       = 26, COBALT     = 27, NICKEL     = 28, COPPER       = 29, ZINC       = 30,
    GALLIUM    = 31, GERMANIUM  = 32, ARSENIC    = 33, SELENIUM     = 34, BROMINE    = 35,
    KRYPTON    = 36, RUBIDIUM   = 37, STRONTIUM  = 38, YTTRIUM      = 39, ZIRCONIUM  = 40,
    NIOBIUM    = 41, MOLYBDENUM = 42, TECHNETIUM = 43, RUTHENIUM    = 44, RHODIUM    = 45,
    PALLADIUM  = 46, SILVER     = 47, CADMIUM    = 48, INDIUM       = 49, TIN        = 50,
    ANTIMONY   = 51, TELLURIUM  = 52, IODINE     = 53, XENON        = 54, CESIUM     = 55,
    BARIUM     = 56, LANTHANUM  = 57, CERIUM     = 58, PRASEODYMIUM = 59, NEODYMIUM  = 60,
    PROMETHIUM = 61, SAMARIUM   = 62, EUROPIUM   = 63, GADOLINIUM   = 64, TERBIUM    = 65,
    DYSPROSIUM = 66, HOLMIUM    = 67, ERBIUM     = 68, THULIUM      = 69, YTTERBIUM  = 70,
    LUTETIUM   = 71, HAFNIUM    = 72, TANTALUM   = 73, TUNGSTEN     = 74, RHENIUM    = 75,
    OSMIUM     = 76, IRIDIUM    = 77, PLATINUM   = 78, GOLD         = 79, MERCURY    = 80,
    THALLIUM   = 81, LEAD       = 82, BISMUTH    = 83, POLONIUM     = 84, ASTATINE   = 85,
    RADON      = 86, FRANCIUM   = 87, RADIUM     = 88, ACTINIUM     = 89, THORIUM    = 90,
    UNK = UINT16_MAX
};

enum class Remoteness : uint8_t
{
    NONE = 0, ALPHA = 1, BETA = 2, GAMMA = 3, DELTA = 4, EPSILON = 5, ZETA = 6, ETA = 7,
    ONE = 11, TWO = 12, THREE = 13, FOUR = 14, FIVE = 15,
    EXTRA = 99,
    UNK = UINT8_MAX
};

enum class Branch : uint8_t
{
    NONE = 0, ONE = 1, TWO = 2, THREE = 3,
    TERMINAL = 254,
    UNK = UINT8_MAX
};

// Ref: https://mmcif.wwpdb.org/dictionaries/mmcif_ma.dic/Items/_struct_conf_type.id.html
enum class Structure : uint8_t
{
    FREE = 0, BEND = 1,
    /* beta-sheet */
    SHEET = 2,
    STRN = 3,
    OTHER = 9,

    /* 10 - 99 for alpha-helix */
    HELX_P       = 10, HELX_OT_P    = 11,
    HELX_LH_27_P = 12, HELX_RH_27_P = 13,
    HELX_LH_3T_P = 14, HELX_RH_3T_P = 15,
    HELX_LH_AL_P = 16, HELX_RH_AL_P = 17,
    HELX_LH_GA_P = 18, HELX_RH_GA_P = 19,
    HELX_LH_OM_P = 20, HELX_RH_OM_P = 21,
    HELX_LH_PI_P = 22, HELX_RH_PI_P = 23,
    HELX_LH_PP_P = 24, HELX_RH_PP_P = 25,
    HELX_LH_OT_P = 26, HELX_RH_OT_P = 27,

    HELX_N       = 50, HELX_OT_N    = 51,
    HELX_LH_N    = 52, HELX_RH_N    = 53,
    HELX_LH_OT_N = 54, HELX_RH_OT_N = 55,
    HELX_LH_A_N  = 56, HELX_RH_A_N  = 57,
    HELX_LH_B_N  = 58, HELX_RH_B_N  = 59,
    HELX_LH_Z_N  = 60, HELX_RH_Z_N  = 61,

    /* 100 - 199 for turn */
    TURN_P      = 100, TURN_OT_P  = 101,
    TURN_TY1P_P = 102, TURN_TY1_P = 103,
    TURN_TY2P_P = 104, TURN_TY2_P = 105,
    TURN_TY3P_P = 106, TURN_TY3_P = 107,
    
    UNK = UINT8_MAX
};

enum class Entity : uint8_t
{
    POLYMER = 0, NONPOLYMER = 1, BRANCHED = 2, MACROLIDE = 3,
    WATER = 4,
    UNK = UINT8_MAX
};
