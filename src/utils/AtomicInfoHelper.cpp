#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "Logger.hpp"

#include <algorithm>
#include <stdexcept>

const std::vector<std::string> AtomicInfoHelper::m_group_atom_class_key_list
{
    "simple_atom_class", "component_atom_class", "structure_atom_class"
};

const std::vector<std::string> AtomicInfoHelper::m_group_bond_class_key_list
{
    "simple_bond_class", "component_bond_class"
};

const std::vector<Residue> AtomicInfoHelper::m_standard_residue_list
{
    Residue::ALA, Residue::ARG, Residue::ASN, Residue::ASP, Residue::CYS,
    Residue::GLN, Residue::GLU, Residue::GLY, Residue::HIS, Residue::ILE,
    Residue::LEU, Residue::LYS, Residue::MET, Residue::PHE, Residue::PRO,
    Residue::SER, Residue::THR, Residue::TRP, Residue::TYR, Residue::VAL
};

const std::vector<Element> AtomicInfoHelper::m_standard_element_list
{
    Element::HYDROGEN, Element::CARBON, Element::NITROGEN, Element::OXYGEN, Element::SULFUR
};

const std::unordered_map<int, Element> AtomicInfoHelper::m_atomic_number_to_element_map
{
    { 1, Element::HYDROGEN},  { 2, Element::HELIUM},      { 3, Element::LITHIUM},
    { 4, Element::BERYLLIUM}, { 5, Element::BORON},       { 6, Element::CARBON},
    { 7, Element::NITROGEN},  { 8, Element::OXYGEN},      { 9, Element::FLUORINE},
    {10, Element::NEON},      {11, Element::SODIUM},      {12, Element::MAGNESIUM},
    {13, Element::ALUMINUM},  {14, Element::SILICON},     {15, Element::PHOSPHORUS},
    {16, Element::SULFUR},    {17, Element::CHLORINE},    {18, Element::ARGON},
    {19, Element::POTASSIUM}, {20, Element::CALCIUM},     {21, Element::SCANDIUM},
    {22, Element::TITANIUM},  {23, Element::VANADIUM},    {24, Element::CHROMIUM},
    {25, Element::MANGANESE}, {26, Element::IRON},        {27, Element::COBALT},
    {28, Element::NICKEL},    {29, Element::COPPER},      {30, Element::ZINC},
    {31, Element::GALLIUM},   {32, Element::GERMANIUM},   {33, Element::ARSENIC},
    {34, Element::SELENIUM},  {35, Element::BROMINE},     {36, Element::KRYPTON},
    {37, Element::RUBIDIUM},  {38, Element::STRONTIUM},   {39, Element::YTTRIUM},
    {40, Element::ZIRCONIUM}, {41, Element::NIOBIUM},     {42, Element::MOLYBDENUM},
    {43, Element::TECHNETIUM},{44, Element::RUTHENIUM},   {45, Element::RHODIUM},
    {46, Element::PALLADIUM}, {47, Element::SILVER},      {48, Element::CADMIUM},
    {49, Element::INDIUM},    {50, Element::TIN},         {51, Element::ANTIMONY},
    {52, Element::TELLURIUM}, {53, Element::IODINE},      {54, Element::XENON},
    {55, Element::CESIUM},    {56, Element::BARIUM},      {57, Element::LANTHANUM},
    {58, Element::CERIUM},    {59, Element::PRASEODYMIUM},{60, Element::NEODYMIUM},
    {61, Element::PROMETHIUM},{62, Element::SAMARIUM},    {63, Element::EUROPIUM},
    {64, Element::GADOLINIUM},{65, Element::TERBIUM},     {66, Element::DYSPROSIUM},
    {67, Element::HOLMIUM},   {68, Element::ERBIUM},      {69, Element::THULIUM},
    {70, Element::YTTERBIUM}, {71, Element::LUTETIUM},    {72, Element::HAFNIUM},
    {73, Element::TANTALUM},  {74, Element::TUNGSTEN},    {75, Element::RHENIUM},
    {76, Element::OSMIUM},    {77, Element::IRIDIUM},     {78, Element::PLATINUM},
    {79, Element::GOLD},      {80, Element::MERCURY},     {81, Element::THALLIUM},
    {82, Element::LEAD},      {83, Element::BISMUTH},     {84, Element::POLONIUM},
    {85, Element::ASTATINE},  {86, Element::RADON},       {87, Element::FRANCIUM},
    {88, Element::RADIUM},    {89, Element::ACTINIUM},    {90, Element::THORIUM}
};

const std::unordered_map<std::string_view, Residue> AtomicInfoHelper::m_residue_map
{
    {"ALA", Residue::ALA}, {"ARG", Residue::ARG}, {"ASN", Residue::ASN}, {"ASP", Residue::ASP},
    {"CYS", Residue::CYS}, {"GLN", Residue::GLN}, {"GLU", Residue::GLU}, {"GLY", Residue::GLY},
    {"HIS", Residue::HIS}, {"ILE", Residue::ILE}, {"LEU", Residue::LEU}, {"LYS", Residue::LYS},
    {"MET", Residue::MET}, {"PHE", Residue::PHE}, {"PRO", Residue::PRO}, {"SER", Residue::SER},
    {"THR", Residue::THR}, {"TRP", Residue::TRP}, {"TYR", Residue::TYR}, {"VAL", Residue::VAL},
    {"A",   Residue::A},   {"C",   Residue::C},   {"G",   Residue::G},   {"U",   Residue::U},
    {"DA",  Residue::DA},  {"DC",  Residue::DC},  {"DG",  Residue::DG},  {"DT",  Residue::DT},
};

const std::unordered_map<std::string_view, Element> AtomicInfoHelper::m_element_map
{
    {"H",  Element::HYDROGEN},  {"HE", Element::HELIUM},      {"LI", Element::LITHIUM},
    {"BE", Element::BERYLLIUM}, {"B",  Element::BORON},       {"C",  Element::CARBON},
    {"N",  Element::NITROGEN},  {"O",  Element::OXYGEN},      {"F",  Element::FLUORINE},
    {"NE", Element::NEON},      {"NA", Element::SODIUM},      {"MG", Element::MAGNESIUM},
    {"AL", Element::ALUMINUM},  {"SI", Element::SILICON},     {"P",  Element::PHOSPHORUS},
    {"S",  Element::SULFUR},    {"CL", Element::CHLORINE},    {"AR", Element::ARGON},
    {"K",  Element::POTASSIUM}, {"CA", Element::CALCIUM},     {"SC", Element::SCANDIUM},
    {"TI", Element::TITANIUM},  {"V",  Element::VANADIUM},    {"CR", Element::CHROMIUM},
    {"MN", Element::MANGANESE}, {"FE", Element::IRON},        {"CO", Element::COBALT},
    {"NI", Element::NICKEL},    {"CU", Element::COPPER},      {"ZN", Element::ZINC},
    {"GA", Element::GALLIUM},   {"GE", Element::GERMANIUM},   {"AS", Element::ARSENIC},
    {"SE", Element::SELENIUM},  {"BR", Element::BROMINE},     {"KR", Element::KRYPTON},
    {"RB", Element::RUBIDIUM},  {"SR", Element::STRONTIUM},   {"Y",  Element::YTTRIUM},
    {"ZR", Element::ZIRCONIUM}, {"NB", Element::NIOBIUM},     {"MO", Element::MOLYBDENUM},
    {"TC", Element::TECHNETIUM},{"RU", Element::RUTHENIUM},   {"RH", Element::RHODIUM},
    {"PD", Element::PALLADIUM}, {"AG", Element::SILVER},      {"CD", Element::CADMIUM},
    {"IN", Element::INDIUM},    {"SN", Element::TIN},         {"SB", Element::ANTIMONY},
    {"TE", Element::TELLURIUM}, {"I",  Element::IODINE},      {"XE", Element::XENON},
    {"CS", Element::CESIUM},    {"BA", Element::BARIUM},      {"LA", Element::LANTHANUM},
    {"CE", Element::CERIUM},    {"PR", Element::PRASEODYMIUM},{"ND", Element::NEODYMIUM},
    {"PM", Element::PROMETHIUM},{"SM", Element::SAMARIUM},    {"EU", Element::EUROPIUM},
    {"GD", Element::GADOLINIUM},{"TB", Element::TERBIUM},     {"DY", Element::DYSPROSIUM},
    {"HO", Element::HOLMIUM},   {"ER", Element::ERBIUM},      {"TM", Element::THULIUM},
    {"YB", Element::YTTERBIUM}, {"LU", Element::LUTETIUM},    {"HF", Element::HAFNIUM},
    {"TA", Element::TANTALUM},  {"W",  Element::TUNGSTEN},    {"RE", Element::RHENIUM},
    {"OS", Element::OSMIUM},    {"IR", Element::IRIDIUM},     {"PT", Element::PLATINUM},
    {"AU", Element::GOLD},      {"HG", Element::MERCURY},     {"TL", Element::THALLIUM},
    {"PB", Element::LEAD},      {"BI", Element::BISMUTH},     {"PO", Element::POLONIUM},
    {"AT", Element::ASTATINE},  {"RN", Element::RADON},       {"FR", Element::FRANCIUM},
    {"RA", Element::RADIUM},    {"AC", Element::ACTINIUM},    {"TH", Element::THORIUM}
};

const std::unordered_map<std::string_view, Spot> AtomicInfoHelper::m_spot_map
{
    {"P",    Spot::P   }, {"OP1",  Spot::OP1 }, {"OP2",  Spot::OP2 }, {"OP3",  Spot::OP3 },
    {"OP1",  Spot::OP1 }, {"OP2",  Spot::OP2 }, {"OP3",  Spot::OP3 },
    {"HOP2", Spot::HOP2}, {"HOP3", Spot::HOP3},
    {"C1'",  Spot::C1p }, {"C2'",  Spot::C2p }, {"C3'",  Spot::C3p }, {"C4'",  Spot::C4p },
    {"C5'",  Spot::C5p },
    {"O2'",  Spot::O2p }, {"O3'",  Spot::O3p }, {"O4'",  Spot::O4p }, {"O5'",  Spot::O5p },
    {"H1'",  Spot::H1p }, {"H2'",  Spot::H2p }, {"H3'",  Spot::H3p }, {"H4'",  Spot::H4p },
    {"H5'",  Spot::H5p }, {"H5''", Spot::H5pp},
    {"HO2'", Spot::HO2p}, {"HO3'", Spot::HO3p},
    {"C2",   Spot::C2  }, {"C4",   Spot::C4  }, {"C5",   Spot::C5  }, {"C6",   Spot::C6  },
    {"C8",   Spot::C8  },
    {"N1",   Spot::N1  }, {"N2",   Spot::N2  }, {"N3",   Spot::N3  }, {"N4",   Spot::N4  },
    {"N6",   Spot::N6  }, {"N7",   Spot::N7  }, {"N9",   Spot::N9  },
    {"O2",   Spot::O2  }, {"O4",   Spot::O4  }, {"O6",   Spot::O6  },
    {"H1",   Spot::H1  }, {"H21",  Spot::H21 }, {"H22",  Spot::H22 }, {"H3",   Spot::H3  },
    {"H41",  Spot::H41 }, {"H42",  Spot::H42 }, {"H5",   Spot::H5  }, {"H6",   Spot::H6  },
    {"H61",  Spot::H61 }, {"H62",  Spot::H62 }, {"H8",   Spot::H8  },

    {"H",    Spot::H   }, {"H2",   Spot::H2  },
    {"HA",   Spot::HA  }, {"HA2",  Spot::HA2 }, {"HA3",  Spot::HA3 },
    {"HB",   Spot::HB  }, {"HB1",  Spot::HB1 }, {"HB2",  Spot::HB2 }, {"HB3",  Spot::HB3 },
    {"HG",   Spot::HG  },
    {"HG1",  Spot::HG1 }, {"HG11", Spot::HG11}, {"HG12", Spot::HG12}, {"HG13", Spot::HG13},
    {"HG2",  Spot::HG2 }, {"HG21", Spot::HG21}, {"HG22", Spot::HG22}, {"HG23", Spot::HG23},
    {"HG3",  Spot::HG3 },
    {"HD1",  Spot::HD1 }, {"HD11", Spot::HD11}, {"HD12", Spot::HD12}, {"HD13", Spot::HD13},
    {"HD2",  Spot::HD2 }, {"HD21", Spot::HD21}, {"HD22", Spot::HD22}, {"HD23", Spot::HD23},
    {"HD3",  Spot::HD3 },
    {"HE",   Spot::HE  }, {"HE1",  Spot::HE1 },
    {"HE2",  Spot::HE2 }, {"HE21", Spot::HE21}, {"HE22", Spot::HE22},
    {"HE3",  Spot::HE3 },
    {"HZ",   Spot::HZ  }, {"HZ1",  Spot::HZ1 }, {"HZ2",  Spot::HZ2 }, {"HZ3",  Spot::HZ3 },
    {"HH",   Spot::HH  }, {"HH11", Spot::HH11}, {"HH12", Spot::HH12},
    {"HH2",  Spot::HH2 }, {"HH21", Spot::HH21}, {"HH22", Spot::HH22},
    {"HXT",  Spot::HXT },

    {"C",    Spot::C   },
    {"CA",   Spot::CA  }, {"CB",   Spot::CB  },
    {"CG",   Spot::CG  }, {"CG1",  Spot::CG1 }, {"CG2",  Spot::CG2 },
    {"CD",   Spot::CD  }, {"CD1",  Spot::CD1 }, {"CD2",  Spot::CD2 },
    {"CE",   Spot::CE  }, {"CE1",  Spot::CE1 }, {"CE2",  Spot::CE2 }, {"CE3",  Spot::CE3 },
    {"CZ",   Spot::CZ  }, {"CZ2",  Spot::CZ2 }, {"CZ3",  Spot::CZ3 },
    {"CH2",  Spot::CH2 },

    {"N",    Spot::N   },
    {"ND1",  Spot::ND1 }, {"ND2",  Spot::ND2 },
    {"NE",   Spot::NE  }, {"NE1",  Spot::NE1 }, {"NE2",  Spot::NE2 },
    {"NZ",   Spot::NZ  },
    {"NH1",  Spot::NH1 }, {"NH2",  Spot::NH2 },

    {"O",    Spot::O   },
    {"OG",   Spot::OG  }, {"OG1",  Spot::OG1 },
    {"OD1",  Spot::OD1 }, {"OD2",  Spot::OD2 },
    {"OE1",  Spot::OE1 }, {"OE2",  Spot::OE2 },
    {"OH",   Spot::OH  },
    {"OXT",  Spot::OXT },

    {"SG",   Spot::SG  }, {"SD",   Spot::SD  }
    
};

const std::unordered_map<std::string_view, Bond> AtomicInfoHelper::m_bond_map
{
    {"OP2_HOP2",  Bond::OP2_HOP2},
    {"OP3_P",     Bond::OP3_P},    {"OP3_HOP3",  Bond::OP3_HOP3},
    {"P_OP1",     Bond::P_OP1},    {"P_OP2",     Bond::P_OP2},    {"P_O5'",     Bond::P_O5p},

    {"C1'_N1",    Bond::C1p_N1},   {"C1'_N9",    Bond::C1p_N9},   {"C1'_H1'",   Bond::C1p_H1p},
    {"C2'_O2'",   Bond::C2p_O2p},  {"C2'_C1'",   Bond::C2p_C1p},  {"C2'_H2'",   Bond::C2p_H2p},
    {"C2'_H2''",  Bond::C2p_H2pp}, {"C3'_O3'",   Bond::C3p_O3p},  {"C3'_C2'",   Bond::C3p_C2p},
    {"C3'_H3'",   Bond::C3p_H3p},  {"C4'_O4'",   Bond::C4p_O4p},  {"C4'_C3'",   Bond::C4p_C3p},
    {"C4'_H4'",   Bond::C4p_H4p},  {"C5'_C4'",   Bond::C5p_C4p},  {"C5'_H5'",   Bond::C5p_H5p},
    {"C5'_H5''",  Bond::C5p_H5pp}, {"O2'_HO2'",  Bond::O2p_HO2p}, {"O3'_HO3'",  Bond::O3p_HO3p},
    {"O4'_C1'",   Bond::O4p_C1p},  {"O5'_C5'",   Bond::O5p_C5p},

    {"C2_H2",     Bond::C2_H2},    {"C2_O2",     Bond::C2_O2},    {"C2_N2",     Bond::C2_N2},
    {"C2_N3",     Bond::C2_N3},
    {"C4_N4",     Bond::C4_N4},    {"C4_O4",     Bond::C4_O4},    {"C4_C5",     Bond::C4_C5},
    {"C5_C4",     Bond::C5_C4},    {"C5_H5",     Bond::C5_H5},    {"C5_C6",     Bond::C5_C6},
    {"C5_C7",     Bond::C5_C7},
    {"C6_N1",     Bond::C6_N1},    {"C6_N6",     Bond::C6_N6},    {"C6_H6",     Bond::C6_H6},
    {"C6_O6",     Bond::C6_O6},
    {"C7_H71",    Bond::C7_H71},   {"C7_H72",    Bond::C7_H72},   {"C7_H73",    Bond::C7_H73},
    {"C8_N7",     Bond::C8_N7},    {"C8_H8",     Bond::C8_H8},    {"N1_H1",     Bond::N1_H1},
    {"N1_C2",     Bond::N1_C2},    {"N1_C6",     Bond::N1_C6},    {"N2_H21",    Bond::N2_H21},
    {"N2_H22",    Bond::N2_H22},   {"N3_C4",     Bond::N3_C4},    {"N3_H3",     Bond::N3_H3},
    {"N4_H41",    Bond::N4_H41},   {"N4_H42",    Bond::N4_H42},   {"N6_H61",    Bond::N6_H61},
    {"N6_H62",    Bond::N6_H62},   {"N7_C5",     Bond::N7_C5},    {"N9_C4",     Bond::N9_C4},
    {"N9_C8",     Bond::N9_C8},

    {"N_CA", Bond::N_CA}, {"N_CD", Bond::N_CD}, {"N_H", Bond::N_H},
    {"N_H2", Bond::N_H2},
    {"CA_C", Bond::CA_C}, {"CA_CB", Bond::CA_CB}, {"CA_HA", Bond::CA_HA},
    {"C_N", Bond::C_N},
    {"C_O", Bond::C_O}, {"C_OXT", Bond::C_OXT},

    {"CB_CG", Bond::CB_CG}, {"CB_OG", Bond::CB_OG}, {"CB_SG", Bond::CB_SG},
    {"CB_CG1", Bond::CB_CG1}, {"CB_CG2", Bond::CB_CG2}, {"CB_OG1", Bond::CB_OG1},
    {"CB_HB", Bond::CB_HB},
    {"CB_HB1", Bond::CB_HB1}, {"CB_HB2", Bond::CB_HB2}, {"CB_HB3", Bond::CB_HB3},
    {"CG_CD", Bond::CG_CD}, {"CG_CD1", Bond::CG_CD1}, {"CG_CD2", Bond::CG_CD2},
    {"CG_HG", Bond::CG_HG}, {"CG_HG2", Bond::CG_HG2}, {"CG_HG3", Bond::CG_HG3},
    {"CG_OD1", Bond::CG_OD1}, {"CG_OD2", Bond::CG_OD2}, {"CG_ND1", Bond::CG_ND1},
    {"CG_ND2", Bond::CG_ND2}, {"CG_SD", Bond::CG_SD},

    {"CG1_CD1", Bond::CG1_CD1}, {"CG1_HG11", Bond::CG1_HG11}, {"CG1_HG12", Bond::CG1_HG12},
    {"CG1_HG13", Bond::CG1_HG13}, {"CG2_HG21", Bond::CG2_HG21}, {"CG2_HG22", Bond::CG2_HG22},
    {"CG2_HG23", Bond::CG2_HG23}, {"CD_CE", Bond::CD_CE}, {"CD_NE", Bond::CD_NE},
    {"CD_NE2", Bond::CD_NE2}, {"CD_OE1", Bond::CD_OE1}, {"CD_OE2", Bond::CD_OE2},
    {"CD_HD2", Bond::CD_HD2}, {"CD_HD3", Bond::CD_HD3}, {"CD1_CE1", Bond::CD1_CE1},
    {"CD1_NE1", Bond::CD1_NE1}, {"CD1_HD1", Bond::CD1_HD1}, {"CD1_HD11", Bond::CD1_HD11},
    {"CD1_HD12", Bond::CD1_HD12}, {"CD1_HD13", Bond::CD1_HD13}, {"CD2_CE2", Bond::CD2_CE2},
    {"CD2_CE3", Bond::CD2_CE3}, {"CD2_NE2", Bond::CD2_NE2}, {"CD2_HD2", Bond::CD2_HD2},
    {"CD2_HD21", Bond::CD2_HD21}, {"CD2_HD22", Bond::CD2_HD22}, {"CD1_HD23", Bond::CD1_HD23},

    {"CE_NZ", Bond::CE_NZ}, {"CE_HE1", Bond::CE_HE1}, {"CE_HE2", Bond::CE_HE2},
    {"CE_HE3", Bond::CE_HE3}, {"CE1_CZ", Bond::CE1_CZ}, {"CE1_NE2", Bond::CE1_NE2},
    {"CE1_HE1", Bond::CE1_HE1}, {"CE2_CZ", Bond::CE2_CZ}, {"CE2_CZ2", Bond::CE2_CZ2},
    {"CE2_HE2", Bond::CE2_HE2}, {"CE3_CZ3", Bond::CE3_CZ3}, {"CE3_HE3", Bond::CE3_HE3},
    {"CZ_HZ", Bond::CZ_HZ}, {"CZ_NH1", Bond::CZ_NH1}, {"CZ_NH2", Bond::CZ_NH2},
    {"CZ_OH", Bond::CZ_OH}, {"CZ2_CH2", Bond::CZ2_CH2}, {"CZ2_HZ2", Bond::CZ2_HZ2},
    {"CZ3_CH2", Bond::CZ3_CH2}, {"CZ3_HZ3", Bond::CZ3_HZ3}, {"CH2_HH2", Bond::CH2_HH2},

    {"ND1_CE1", Bond::ND1_CE1}, {"ND1_HD1", Bond::ND1_HD1},
    {"ND2_HD21", Bond::ND2_HD21}, {"ND2_HD22", Bond::ND2_HD22},
    {"NE_CZ", Bond::NE_CZ}, {"NE_HE", Bond::NE_HE}, {"NE1_CE2", Bond::NE1_CE2},
    {"NE1_HE1", Bond::NE1_HE1}, {"NE2_HE2", Bond::NE2_HE2}, {"NE2_HE21", Bond::NE2_HE21},
    {"NE2_HE22", Bond::NE2_HE22}, {"NZ_HZ1", Bond::NZ_HZ1}, {"NZ_HZ2", Bond::NZ_HZ2},
    {"NZ_HZ3", Bond::NZ_HZ3}, {"NH1_HH11", Bond::NH1_HH11}, {"NH1_HH12", Bond::NH1_HH12},
    {"NH2_HH21", Bond::NH2_HH21}, {"NH2_HH22", Bond::NH2_HH22},

    {"O_H1", Bond::O_H1}, {"O_H2", Bond::O_H2}, {"OG_HG", Bond::OG_HG},
    {"OG1_HG1", Bond::OG1_HG1}, {"OD2_HD2", Bond::OD2_HD2}, {"OE2_HE2", Bond::OE2_HE2},
    {"OH_HH", Bond::OH_HH}, {"OXT_HXT", Bond::OXT_HXT}, {"SD_CE", Bond::SD_CE}
};

const std::unordered_map<std::string_view, Structure> AtomicInfoHelper::m_structure_map
{
    {" ", Structure::FREE}, {"BEND", Structure::BEND},
    {"STRN", Structure::STRN},
    {"OTHER", Structure::OTHER},
    {"HELX_P", Structure::HELX_P},
    {"TURN_P", Structure::TURN_P}
};

const std::unordered_map<std::string_view, Entity> AtomicInfoHelper::m_entity_map
{
    {"polymer",  Entity::POLYMER},  {"non-polymer", Entity::NONPOLYMER},
    {"branched", Entity::BRANCHED}, {"macrolide",   Entity::MACROLIDE},
    {"water",    Entity::WATER}
};

const std::unordered_map<Residue, std::string> AtomicInfoHelper::m_residue_label_map
{
    {Residue::ALA, "ALA"}, {Residue::ARG, "ARG"}, {Residue::ASN, "ASN"}, {Residue::ASP, "ASP"},
    {Residue::CYS, "CYS"}, {Residue::GLN, "GLN"}, {Residue::GLU, "GLU"}, {Residue::GLY, "GLY"},
    {Residue::HIS, "HIS"}, {Residue::ILE, "ILE"}, {Residue::LEU, "LEU"}, {Residue::LYS, "LYS"},
    {Residue::MET, "MET"}, {Residue::PHE, "PHE"}, {Residue::PRO, "PRO"}, {Residue::SER, "SER"},
    {Residue::THR, "THR"}, {Residue::TRP, "TRP"}, {Residue::TYR, "TYR"}, {Residue::VAL, "VAL"},
    {Residue::A,   "A"  }, {Residue::C,   "C"  }, {Residue::G,   "G"  }, {Residue::U,   "U"  },
    {Residue::DA,  "DA" }, {Residue::DC,  "DC" }, {Residue::DG,  "DG" }, {Residue::DT,  "DT" },
    {Residue::UNK, "UNK"}
};

const std::unordered_map<Element, std::string> AtomicInfoHelper::m_element_label_map
{
    {Element::HYDROGEN,   "H" }, {Element::HELIUM,      "He"}, {Element::LITHIUM,    "Li"},
    {Element::BERYLLIUM,  "Be"}, {Element::BORON,       "B" }, {Element::CARBON,     "C" },
    {Element::NITROGEN,   "N" }, {Element::OXYGEN,      "O" }, {Element::FLUORINE,   "F" },
    {Element::NEON,       "Ne"}, {Element::SODIUM,      "Na"}, {Element::MAGNESIUM,  "Mg"},
    {Element::ALUMINUM,   "Al"}, {Element::SILICON,     "Si"}, {Element::PHOSPHORUS, "P" },
    {Element::SULFUR,     "S" }, {Element::CHLORINE,    "Cl"}, {Element::ARGON,      "Ar"},
    {Element::POTASSIUM,  "K" }, {Element::CALCIUM,     "Ca"}, {Element::SCANDIUM,   "Sc"},
    {Element::TITANIUM,   "Ti"}, {Element::VANADIUM,    "V" }, {Element::CHROMIUM,   "Cr"},
    {Element::MANGANESE,  "Mn"}, {Element::IRON,        "Fe"}, {Element::COBALT,     "Co"},
    {Element::NICKEL,     "Ni"}, {Element::COPPER,      "Cu"}, {Element::ZINC,       "Zn"},
    {Element::GALLIUM,    "Ga"}, {Element::GERMANIUM,   "Ge"}, {Element::ARSENIC,    "As"},
    {Element::SELENIUM,   "Se"}, {Element::BROMINE,     "Br"}, {Element::KRYPTON,    "Kr"},
    {Element::RUBIDIUM,   "Rb"}, {Element::STRONTIUM,   "Sr"}, {Element::YTTRIUM,    "Y" },
    {Element::ZIRCONIUM,  "Zr"}, {Element::NIOBIUM,     "Nb"}, {Element::MOLYBDENUM, "Mo"},
    {Element::TECHNETIUM, "Tc"}, {Element::RUTHENIUM,   "Ru"}, {Element::RHODIUM,    "Rh"},
    {Element::PALLADIUM,  "Pd"}, {Element::SILVER,      "Ag"}, {Element::CADMIUM,    "Cd"},
    {Element::INDIUM,     "In"}, {Element::TIN,         "Sn"}, {Element::ANTIMONY,   "Sb"},
    {Element::TELLURIUM,  "Te"}, {Element::IODINE,      "I" }, {Element::XENON,      "Xe"},
    {Element::CESIUM,     "Cs"}, {Element::BARIUM,      "Ba"}, {Element::LANTHANUM,  "La"},
    {Element::CERIUM,     "Ce"}, {Element::PRASEODYMIUM,"Pr"}, {Element::NEODYMIUM,  "Nd"},
    {Element::PROMETHIUM, "Pm"}, {Element::SAMARIUM,    "Sm"}, {Element::EUROPIUM,   "Eu"},
    {Element::GADOLINIUM, "Gd"}, {Element::TERBIUM,     "Tb"}, {Element::DYSPROSIUM, "Dy"},
    {Element::HOLMIUM,    "Ho"}, {Element::ERBIUM,      "Er"}, {Element::THULIUM,    "Tm"},
    {Element::YTTERBIUM,  "Yb"}, {Element::LUTETIUM,    "Lu"}, {Element::HAFNIUM,    "Hf"},
    {Element::TANTALUM,   "Ta"}, {Element::TUNGSTEN,    "W" }, {Element::RHENIUM,    "Re"},
    {Element::OSMIUM,     "Os"}, {Element::IRIDIUM,     "Ir"}, {Element::PLATINUM,   "Pt"},
    {Element::GOLD,       "Au"}, {Element::MERCURY,     "Hg"}, {Element::THALLIUM,   "Tl"},
    {Element::LEAD,       "Pb"}, {Element::BISMUTH,     "Bi"}, {Element::POLONIUM,   "Po"},
    {Element::ASTATINE,   "At"}, {Element::RADON,       "Rn"}, {Element::FRANCIUM,   "Fr"},
    {Element::RADIUM,     "Ra"}, {Element::ACTINIUM,    "Ac"}, {Element::THORIUM,    "Th"},
    {Element::UNK,       "UNK"}
};

const std::unordered_map<Element, int> AtomicInfoHelper::m_element_color_map
{
    {Element::HYDROGEN, 921}, {Element::CARBON,   633}, {Element::NITROGEN,  418},
    {Element::OXYGEN,   862}, {Element::PHOSPHORUS, 2},
    {Element::SULFUR,   619}, {Element::CALCIUM,    1},
    {Element::ZINC,       1}, {Element::SODIUM,     1}, {Element::MAGNESIUM,   1},
    {Element::IRON,       1}, {Element::CHLORINE,   1}, {Element::UNK,         1}
};

const std::unordered_map<Residue, int> AtomicInfoHelper::m_residue_color_map
{
    {Residue::ALA, 632}, {Residue::ARG, 633}, {Residue::ASN, 634}, {Residue::ASP, 635},
    {Residue::CYS, 800}, {Residue::GLN, 801}, {Residue::GLU, 802}, {Residue::GLY, 803},
    {Residue::HIS, 416}, {Residue::ILE, 417}, {Residue::LEU, 418}, {Residue::LYS, 419},
    {Residue::MET, 600}, {Residue::PHE, 601}, {Residue::PRO, 602}, {Residue::SER, 603},
    {Residue::THR, 880}, {Residue::TRP, 881}, {Residue::TYR, 882}, {Residue::VAL, 883},

    {Residue::A,     2}, {Residue::C,     3}, {Residue::G,     4}, {Residue::U,     5},
    {Residue::DA,    2}, {Residue::DC,    3}, {Residue::DG,    4}, {Residue::DT,    5}
};

const std::unordered_map<Element, int> AtomicInfoHelper::m_element_marker_map
{
    {Element::HYDROGEN,   5}, {Element::HELIUM,    2}, {Element::LITHIUM,     3},
    {Element::BERYLLIUM, 30}, {Element::BORON,    28}, {Element::CARBON,     53},
    {Element::NITROGEN,  55}, {Element::OXYGEN,   59}, {Element::FLUORINE,   35},
    {Element::NEON,      36}, {Element::SODIUM,   29}, {Element::MAGNESIUM,  22},
    {Element::ALUMINUM,  34}, {Element::SILICON,  41}, {Element::PHOSPHORUS, 21},
    {Element::SULFUR,    33}, {Element::CHLORINE, 39}, {Element::ARGON,      20},
    {Element::POTASSIUM, 31}, {Element::CALCIUM,  45}, {Element::SCANDIUM,   36},
    {Element::TITANIUM,  37}, {Element::VANADIUM, 38}, {Element::CHROMIUM,   40},
    {Element::MANGANESE, 43}, {Element::IRON,     47}, {Element::COBALT,     44},
    {Element::NICKEL,    46}, {Element::COPPER,   49}, {Element::ZINC,       48},
    {Element::CADMIUM,   49},
    {Element::UNK,        1}
};

const std::unordered_map<Residue, int> AtomicInfoHelper::m_residue_marker_map
{
    {Residue::ALA, 20}, {Residue::ARG, 21}, {Residue::ASN, 22}, {Residue::ASP, 23},
    {Residue::CYS, 24}, {Residue::GLN, 25}, {Residue::GLU, 26}, {Residue::GLY, 27},
    {Residue::HIS, 28}, {Residue::ILE, 29}, {Residue::LEU, 30}, {Residue::LYS, 31},
    {Residue::MET, 32}, {Residue::PHE, 33}, {Residue::PRO, 34}, {Residue::SER, 35},
    {Residue::THR, 36}, {Residue::TRP, 37}, {Residue::TYR, 38}, {Residue::VAL, 39},

    {Residue::A,    2}, {Residue::C,    3}, {Residue::G,    4}, {Residue::U,    5},
    {Residue::DA,   2}, {Residue::DC,   3}, {Residue::DG,   4}, {Residue::DT,   5}
};

int AtomicInfoHelper::GetAtomicNumber(Element element)
{
    return static_cast<int>(element);
}

size_t AtomicInfoHelper::GetGroupAtomClassCount(void)
{
    return m_group_atom_class_key_list.size();
}

size_t AtomicInfoHelper::GetGroupBondClassCount(void)
{
    return m_group_bond_class_key_list.size();
}

size_t AtomicInfoHelper::GetStandardResidueCount(void)
{
    return m_standard_residue_list.size();
}

size_t AtomicInfoHelper::GetElementCount(void)
{
    return m_element_map.size();
}

const std::string & AtomicInfoHelper::GetGroupAtomClassKey(size_t class_id)
{
    if (class_id >= m_group_atom_class_key_list.size())
    {
        throw std::out_of_range("class_id out of range");
    }
    return m_group_atom_class_key_list.at(class_id);
}

const std::string & AtomicInfoHelper::GetGroupBondClassKey(size_t class_id)
{
    if (class_id >= m_group_bond_class_key_list.size())
    {
        throw std::out_of_range("class_id out of range");
    }
    return m_group_bond_class_key_list.at(class_id);
}

const std::string & AtomicInfoHelper::GetSimpleAtomClassKey(void)
{
    return m_group_atom_class_key_list.at(0);
}

const std::string & AtomicInfoHelper::GetComponentAtomClassKey(void)
{
    return m_group_atom_class_key_list.at(1);
}

const std::string & AtomicInfoHelper::GetStructureAtomClassKey(void)
{
    return m_group_atom_class_key_list.at(2);
}

const std::string & AtomicInfoHelper::GetSimpleBondClassKey(void)
{
    return m_group_bond_class_key_list.at(0);
}

const std::string & AtomicInfoHelper::GetComponentBondClassKey(void)
{
    return m_group_bond_class_key_list.at(1);
}

const std::vector<Residue> & AtomicInfoHelper::GetStandardResidueList(void)
{
    return m_standard_residue_list;
}

const std::vector<Element> & AtomicInfoHelper::GetStandardElementList(void)
{
    return m_standard_element_list;
}

const std::unordered_map<std::string_view, Residue> & AtomicInfoHelper::GetResidueMap(void)
{
    return m_residue_map;
}

const std::unordered_map<std::string_view, Element> & AtomicInfoHelper::GetElementMap(void)
{
    return m_element_map;
}

const std::unordered_map<std::string_view, Spot> & AtomicInfoHelper::GetSpotMap(void)
{
    return m_spot_map;
}

const std::unordered_map<std::string_view, Bond> & AtomicInfoHelper::GetBondMap(void)
{
    return m_bond_map;
}

const std::unordered_map<Element, std::string> & AtomicInfoHelper::GetElementLabelMap(void)
{
    return m_element_label_map;
}

bool AtomicInfoHelper::IsStandardElement(Element element)
{
    return std::find(m_standard_element_list.begin(), m_standard_element_list.end(), element)
            != m_standard_element_list.end();
}

bool AtomicInfoHelper::IsStandardResidue(Residue residue)
{
    return std::find(m_standard_residue_list.begin(), m_standard_residue_list.end(), residue)
            != m_standard_residue_list.end();
}

Residue AtomicInfoHelper::GetResidueFromString(const std::string & name, bool verbose)
{
    thread_local static std::unordered_map<std::string, int> unknown_name_count_list;
    if (m_residue_map.find(name) == m_residue_map.end())
    {
        if (unknown_name_count_list.find(name) == unknown_name_count_list.end())
        {
            if (verbose)
            {
                Logger::Log(LogLevel::Warning, 
                    "AtomicInfoHelper::GetResidueFromString - Unknown string: " + name);
            }
            unknown_name_count_list[name] = 1;
        }
        else
        {
            unknown_name_count_list[name]++;
        }
        return Residue::UNK;
    }
    return m_residue_map.at(name);
}

Element AtomicInfoHelper::GetElementFromString(const std::string & name)
{
    thread_local static std::unordered_map<std::string, int> unknown_name_count_list;
    if (m_element_map.find(name) == m_element_map.end())
    {
        if (unknown_name_count_list.find(name) == unknown_name_count_list.end())
        {
            Logger::Log(LogLevel::Warning, 
                        "AtomicInfoHelper::GetElementFromString - Unknown string: " + name);
            unknown_name_count_list[name] = 1;
        }
        else
        {
            unknown_name_count_list[name]++;
        }
        return Element::UNK;
    }
    return m_element_map.at(name);
}

Element AtomicInfoHelper::GetElementFromAtomicNumber(int atomic_number)
{
    thread_local static std::unordered_map<int, int> unknown_atomic_number_count_list;
    if (m_atomic_number_to_element_map.find(atomic_number) == m_atomic_number_to_element_map.end())
    {
        if (unknown_atomic_number_count_list.find(atomic_number) == unknown_atomic_number_count_list.end())
        {
            Logger::Log(LogLevel::Warning, 
                "AtomicInfoHelper::GetElementFromAtomicNumber - Unknown atomic number: "
                + std::to_string(atomic_number));
                unknown_atomic_number_count_list[atomic_number] = 1;
        }
        else
        {
            unknown_atomic_number_count_list[atomic_number]++;
        }
        return Element::UNK;
    }
    return m_atomic_number_to_element_map.at(atomic_number);
}

Spot AtomicInfoHelper::GetSpotFromString(const std::string & name, bool verbose)
{
    thread_local static std::unordered_map<std::string, int> unknown_name_count_list;
    if (m_spot_map.find(name) == m_spot_map.end())
    {
        if (unknown_name_count_list.find(name) == unknown_name_count_list.end())
        {
            if (verbose)
            {
                Logger::Log(LogLevel::Warning, 
                    "AtomicInfoHelper::GetSpotFromString - Unknown string: " + name);
            }
            unknown_name_count_list[name] = 1;
        }
        else
        {
            unknown_name_count_list[name]++;
        }
        return Spot::UNK;
    }
    return m_spot_map.at(name);
}

Bond AtomicInfoHelper::GetBondFromString(const std::string & name, bool verbose)
{
    thread_local static std::unordered_map<std::string, int> unknown_name_count_list;
    if (m_bond_map.find(name) == m_bond_map.end())
    {
        if (unknown_name_count_list.find(name) == unknown_name_count_list.end())
        {
            if (verbose)
            {
                Logger::Log(LogLevel::Warning, 
                    "AtomicInfoHelper::GetBondFromString - Unknown string: " + name);
            }
            unknown_name_count_list[name] = 1;
        }
        else
        {
            unknown_name_count_list[name]++;
        }
        return Bond::UNK;
    }
    return m_bond_map.at(name);
}

Structure AtomicInfoHelper::GetStructureFromString(const std::string & name)
{
    thread_local static std::unordered_map<std::string, int> unknown_name_count_list;
    if (m_structure_map.find(name) == m_structure_map.end())
    {
        if (unknown_name_count_list.find(name) == unknown_name_count_list.end())
        {
            Logger::Log(LogLevel::Warning, 
                        "AtomicInfoHelper::GetStructureFromString - Unknown string: " + name);
            unknown_name_count_list[name] = 1;
        }
        else
        {
            unknown_name_count_list[name]++;
        }
        return Structure::UNK;
    }
    return m_structure_map.at(name);
}

Entity AtomicInfoHelper::GetEntityFromString(const std::string & name)
{
    thread_local static std::unordered_map<std::string, int> unknown_name_count_list;
    if (m_entity_map.find(name) == m_entity_map.end())
    {
        if (unknown_name_count_list.find(name) == unknown_name_count_list.end())
        {
            Logger::Log(LogLevel::Warning, 
                "AtomicInfoHelper::GetEntityFromString - Unknown string: " + name);
            unknown_name_count_list[name] = 1;
        }
        else
        {
            unknown_name_count_list[name]++;
        }
        return Entity::UNK;
    }
    return m_entity_map.at(name);
}

const std::string & AtomicInfoHelper::GetLabel(Residue residue)
{
    static std::string unk_label{"?"};
    if (m_residue_label_map.find(residue) == m_residue_label_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "AtomicInfoHelper::GetLabel - Unknown Residue: "
            + std::to_string(static_cast<int>(residue)));
        return unk_label;
    }
    return m_residue_label_map.at(residue);
}

const std::string & AtomicInfoHelper::GetLabel(Element element)
{
    static std::string unk_label{"?"};
    if (m_element_label_map.find(element) == m_element_label_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "AtomicInfoHelper::GetLabel - Unknown Element: "
            + std::to_string(static_cast<int>(element)));
        return unk_label;
    }
    return m_element_label_map.at(element);
}

short AtomicInfoHelper::GetDisplayColor(Element element)
{
    if (m_element_color_map.find(element) == m_element_color_map.end()) return 1;
    return static_cast<short>(m_element_color_map.at(element));
}

short AtomicInfoHelper::GetDisplayColor(Residue residue)
{
    if (m_residue_color_map.find(residue) == m_residue_color_map.end()) return 1;
    return static_cast<short>(m_residue_color_map.at(residue));
}

short AtomicInfoHelper::GetDisplayMarker(Element element)
{
    if (m_element_marker_map.find(element) == m_element_marker_map.end()) return 121;
    return static_cast<short>(m_element_marker_map.at(element));
}

short AtomicInfoHelper::GetDisplayMarker(Residue residue)
{
    if (m_residue_marker_map.find(residue) == m_residue_marker_map.end()) return 1;
    return static_cast<short>(m_residue_marker_map.at(residue));
}
