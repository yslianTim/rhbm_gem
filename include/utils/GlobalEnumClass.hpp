#pragma once

#include <cstdint>

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
    HYDROGEN = 1, CARBON = 6, NITROGEN = 7, OXYGEN = 8,
    SODIUM = 11, MAGNESIUM = 12, SULFUR = 16, CHLORINE = 17,
    CALCIUM = 20, IRON = 26, ZINC = 30,
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
    FREE = 0, BEND = 1, STRN = 2, OTHER = 3,
    HELX_P = 10,
    TURN_P = 100,
    UNK = UINT8_MAX
};