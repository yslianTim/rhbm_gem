#pragma once

#include <cstdint>

enum class Residue : uint8_t
{
    UNK =  0,

    ALA =  1, ARG =  2, ASN =  3, ASP =  4, CYS =  5,
    GLN =  6, GLU =  7, GLY =  8, HIS =  9, ILE = 10,
    LEU = 11, LYS = 12, MET = 13, PHE = 14, PRO = 15,
    SER = 16, THR = 17, TRP = 18, TYR = 19, VAL = 20,
    CSX = 21,

    A   = 25, G   = 26, C   = 27, U   = 28, T   = 29

};

enum class Spot : uint32_t
{
    UNK = 0,

    // Phosphoric acid
    P = 1,
    OP1 = 2, OP2 = 3, OP3 = 4,
    HOP2 = 5, HOP3 = 6,

    // Pentose sugar (Ribose, Deoxyribose)
    C1p  = 11, C2p  = 12, C3p = 13, C4p = 14, C5p = 15,
    O2p  = 22, O3p  = 23, O4p = 24, O5p = 25,
    H1p  = 31, H2p  = 32, H3p = 33, H4p = 32, H5p = 35, H5pp = 36,
    HO2p = 37, HO3p = 38,

    // Nucleotide base
    C2  = 42, C4  = 44, C5  = 45, C6 = 46, C8 = 48,
    N1  = 51, N2  = 52, N3  = 53, N4 = 54, N6 = 56, N7 = 57, N9 = 59,
    O2  = 62, O4  = 64, O6  = 66,
    H1  = 71, H21 = 72, H22 = 73, H3 = 74,  // there's H2 but shared with reisude's one
    H41 = 75, H42 = 76, H5  = 77,
    H6  = 78, H61 = 79, H62 = 80, H8 = 81,

    // Amino Acid components
    H    = 100, H2   = 120,
    HA   = 200, HA2  = 220, HA3  = 230,
    HB   = 300, HB1  = 310, HB2  = 320, HB3  = 330,
    HG   = 400, HG1  = 410, HG11 = 411, HG12 = 412, HG13 = 413,
    HG2  = 420, HG21 = 421, HG22 = 422, HG23 = 423, HG3  = 430,
    HD1  = 510, HD11 = 511, HD12 = 512, HD13 = 513,
    HD2  = 520, HD21 = 521, HD22 = 522, HD23 = 523, HD3  = 530,
    HE   = 600, HE1  = 610, HE2  = 620, HE21 = 621, HE22 = 622, HE3  = 630,
    HZ   = 700, HZ1  = 710, HZ2  = 720, HZ3  = 730,
    HH   = 800, HH11 = 811, HH12 = 812, HH2  = 820, HH21 = 821, HH22 = 822,
    HXT  = 999,

    C   = 6100,
    CA  = 6200,
    CB  = 6300,
    CG  = 6400, CG1 = 6410, CG2 = 6420,
    CD  = 6500, CD1 = 6510, CD2 = 6520,
    CE  = 6600, CE1 = 6610, CE2 = 6620, CE3 = 6630,
    CZ  = 6700, CZ2 = 6720, CZ3 = 6730,
    CH2 = 6820,

    N   = 7100,
    ND1 = 7510, ND2 = 7520,
    NE  = 7600, NE1 = 7610, NE2 = 7620,
    NZ  = 7700, NH1 = 7710, NH2 = 7720,

    O   = 8100,
    OG  = 8400, OG1 = 8410,
    OD1 = 8510, OD2 = 8520,
    OE1 = 8610, OE2 = 8620,
    OH  = 8800,
    OXT = 8999,

    SG =  9400,
    SD =  9500

};

enum class Element : uint8_t
{
    UNK = 0,
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
    RADON      = 86, FRANCIUM   = 87, RADIUM     = 88, ACTINIUM     = 89, THORIUM    = 90
};

enum class Remoteness : uint8_t
{
    UNK = 0,
    NONE = 1, ALPHA = 2, BETA = 3, GAMMA = 4, DELTA = 5, EPSILON = 6, ZETA = 7, ETA = 8,
    EXTRA = 9,
    ACID = 10,         // phosphoric acid
    PENTOSE = 11,      // pentose sugar (Ribose, Deoxyribose)
    BASE = 12          // nucleotide base
};

enum class Branch : uint8_t
{
    UNK = 0,
    NONE = 1, ONE = 2, TWO = 3, THREE = 4, FOUR = 5,
    FIVE = 6, SIX = 7, SEVEN = 8, EIGHT = 9, NINE = 10,
    TERMINAL = 11
};

// Ref: https://mmcif.wwpdb.org/dictionaries/mmcif_ma.dic/Items/_struct_conf_type.id.html
enum class Structure : uint8_t
{
    UNK = 0,
    FREE = 1, BEND = 2,
    /* beta-sheet */
    SHEET = 3,
    STRN = 4,
    OTHER = 5,

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
    TURN_TY3P_P = 106, TURN_TY3_P = 107
};

enum class Entity : uint8_t
{
    UNK = 0,
    POLYMER = 1, NONPOLYMER = 2, BRANCHED = 3, MACROLIDE = 4,
    WATER = 5
};
