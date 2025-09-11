#pragma once

#include <cstdint>

using GroupKey = uint64_t;

enum class Residue : uint8_t
{
    UNK =  0,

    ALA =  1, ARG =  2, ASN =  3, ASP =  4, CYS =  5,
    GLN =  6, GLU =  7, GLY =  8, HIS =  9, ILE = 10,
    LEU = 11, LYS = 12, MET = 13, PHE = 14, PRO = 15,
    SER = 16, THR = 17, TRP = 18, TYR = 19, VAL = 20,
    A   = 21, C   = 22, G   = 23, U   = 24,
    DA  = 25, DC  = 26, DG  = 27, DT  = 28
};

enum class Spot : uint16_t
{
    UNK = 0,

    // Phosphoric acid
    P = 1,
    OP1 = 2, OP2 = 3, OP3 = 4,
    HOP2 = 5, HOP3 = 6,

    // Pentose sugar (Ribose, Deoxyribose)
    C1p  = 11, C2p  = 12, C3p = 13, C4p = 14, C5p = 15,
    O2p  = 22, O3p  = 23, O4p = 24, O5p = 25,
    H1p  = 31, H2p  = 32, H3p = 33, H4p = 32, H5p = 35,
    H2pp = 36, H5pp = 37,
    HO2p = 38, HO3p = 39,

    // Nucleotide base
    C2  = 42, C4  = 44, C5  = 45, C6 = 46, C7 = 47, C8 = 48,
    N1  = 51, N2  = 52, N3  = 53, N4 = 54, N6 = 56, N7 = 57, N9 = 59,
    O2  = 62, O4  = 64, O6  = 66,
    H1  = 71, H21 = 72, H22 = 73, H3 = 74,  // there's H2 but shared with reisude's one
    H41 = 75, H42 = 76, H5  = 77,
    H6  = 78, H61 = 79, H62 = 80, H8 = 81,
    H71 = 82, H72 = 83, H73 = 84,

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

enum class Bond : uint16_t
{
    UNK = 0,

    // Phosphoric acid
    OP2_HOP2,
    OP3_P, OP3_HOP3,
    P_OP1, P_OP2, P_O5p,

    // Pentose sugar (Ribose, Deoxyribose)
    C1p_N1, C1p_N9, C1p_H1p,
    C2p_O2p, C2p_C1p, C2p_H2p, C2p_H2pp,
    C3p_O3p, C3p_C2p, C3p_H3p,
    C4p_O4p, C4p_C3p, C4p_H4p,
    C5p_C4p, C5p_H5p, C5p_H5pp,
    O2p_HO2p, O3p_HO3p, O4p_C1p, O5p_C5p,
    
    // Nucleotide base
    C2_H2, C2_O2, C2_N2, C2_N3,
    C4_N4, C4_O4, C4_C5,
    C5_C4, C5_H5, C5_C6, C5_C7,
    C6_N1, C6_N6, C6_H6, C6_O6,
    C7_H71, C7_H72, C7_H73,
    C8_N7, C8_H8,
    N1_H1, N1_C2, N1_C6,
    N2_H21, N2_H22,
    N3_C4, N3_H3,
    N4_H41, N4_H42,
    N6_H61, N6_H62,
    N7_C5,
    N9_C4, N9_C8,

    
    // Amino Acid components bond
    N_CA, N_CD, N_H, N_H2,
    CA_C, CA_CB, CA_HA,
    C_N,   // peptide connection 
    C_O, C_OXT,

    CB_CG, CB_OG, CB_SG, CB_HB1, CB_HB2, CB_HB3,
    CG_CD, CG_CD1, CG_CD2, CG_HG, CG_HG2, CG_HG3, CG_OD1, CG_OD2, CG_ND1, CG_ND2, CG_SD,
    CG1_CD1, CG1_HG11, CG1_HG12, CG1_HG13,
    CG2_HG21, CG2_HG22, CG2_HG23,
    CD_CE, CD_NE, CD_NE2, CD_OE1, CD_OE2, CD_HD2, CD_HD3,
    CD1_CE1, CD1_NE1, CD1_HD1, CD1_HD11, CD1_HD12, CD1_HD13,
    CD2_CE2, CD2_CE3, CD2_NE2, CD2_HD2, CD2_HD21, CD2_HD22, CD1_HD23,
    CE_NZ, CE_HE1, CE_HE2, CE_HE3,
    CE1_CZ, CE1_NE2, CE1_HE1,
    CE2_CZ, CE2_CZ2, CE2_HE2,
    CE3_CZ3, CE3_HE3,
    CZ_HZ, CZ_NH1, CZ_NH2, CZ_OH,
    CZ2_CH2, CZ2_HZ2,
    CZ3_CH2, CZ3_HZ3,
    CH2_HH2,
    
    ND1_CE1, ND1_HD1, ND2_HD21, ND2_HD22,
    NE_CZ, NE_HE,
    NE1_CE2, NE1_HE1,
    NE2_HE2, NE2_HE21, NE2_HE22,
    NZ_HZ1, NZ_HZ2, NZ_HZ3,
    NH1_HH11, NH1_HH12, NH2_HH21, NH2_HH22,

    O_H1, O_H2,
    OG_HG,
    OG1_HG1,
    OD2_HD2,
    OE2_HE2,
    OH_HH,
    OXT_HXT,

    SD_CE

};

// Element value defined as atomic number
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
