#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "Logger.hpp"

#include <algorithm>
#include <stdexcept>

const std::vector<std::string> AtomicInfoHelper::m_group_class_key_list
{
    "element_class", "residue_class", "structure_class"
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

const std::vector<Remoteness> AtomicInfoHelper::m_standard_remoteness_list
{
    Remoteness::NONE, Remoteness::ALPHA, Remoteness::BETA, Remoteness::GAMMA,
    Remoteness::DELTA, Remoteness::EPSILON, Remoteness::ZETA, Remoteness::ETA
};

const std::vector<Branch> AtomicInfoHelper::m_standard_branch_list
{
    Branch::NONE, Branch::ONE, Branch::TWO, Branch::THREE
};

const std::unordered_map<Element, int> AtomicInfoHelper::m_atomic_number_map
{
    {Element::HYDROGEN,   1}, {Element::HELIUM,       2}, {Element::LITHIUM,     3},
    {Element::BERYLLIUM,  4}, {Element::BORON,        5}, {Element::CARBON,      6},
    {Element::NITROGEN,   7}, {Element::OXYGEN,       8}, {Element::FLUORINE,    9},
    {Element::NEON,      10}, {Element::SODIUM,      11}, {Element::MAGNESIUM,  12},
    {Element::ALUMINUM,  13}, {Element::SILICON,     14}, {Element::PHOSPHORUS, 15},
    {Element::SULFUR,    16}, {Element::CHLORINE,    17}, {Element::ARGON,      18},
    {Element::POTASSIUM, 19}, {Element::CALCIUM,     20}, {Element::SCANDIUM,   21},
    {Element::TITANIUM,  22}, {Element::VANADIUM,    23}, {Element::CHROMIUM,   24},
    {Element::MANGANESE, 25}, {Element::IRON,        26}, {Element::COBALT,     27},
    {Element::NICKEL,    28}, {Element::COPPER,      29}, {Element::ZINC,       30},
    {Element::GALLIUM,   31}, {Element::GERMANIUM,   32}, {Element::ARSENIC,    33},
    {Element::SELENIUM,  34}, {Element::BROMINE,     35}, {Element::KRYPTON,    36},
    {Element::RUBIDIUM,  37}, {Element::STRONTIUM,   38}, {Element::YTTRIUM,    39},
    {Element::ZIRCONIUM, 40}, {Element::NIOBIUM,     41}, {Element::MOLYBDENUM, 42},
    {Element::TECHNETIUM,43}, {Element::RUTHENIUM,   44}, {Element::RHODIUM,    45},
    {Element::PALLADIUM, 46}, {Element::SILVER,      47}, {Element::CADMIUM,    48},
    {Element::INDIUM,    49}, {Element::TIN,         50}, {Element::ANTIMONY,   51},
    {Element::TELLURIUM, 52}, {Element::IODINE,      53}, {Element::XENON,      54},
    {Element::CESIUM,    55}, {Element::BARIUM,      56}, {Element::LANTHANUM,  57},
    {Element::CERIUM,    58}, {Element::PRASEODYMIUM,59}, {Element::NEODYMIUM,  60},
    {Element::PROMETHIUM,61}, {Element::SAMARIUM,    62}, {Element::EUROPIUM,   63},
    {Element::GADOLINIUM,64}, {Element::TERBIUM,     65}, {Element::DYSPROSIUM, 66},
    {Element::HOLMIUM,   67}, {Element::ERBIUM,      68}, {Element::THULIUM,    69},
    {Element::YTTERBIUM, 70}, {Element::LUTETIUM,    71}, {Element::HAFNIUM,    72},
    {Element::TANTALUM,  73}, {Element::TUNGSTEN,    74}, {Element::RHENIUM,    75},
    {Element::OSMIUM,    76}, {Element::IRIDIUM,     77}, {Element::PLATINUM,   78},
    {Element::GOLD,      79}, {Element::MERCURY,     80}, {Element::THALLIUM,   81},
    {Element::LEAD,      82}, {Element::BISMUTH,     83}, {Element::POLONIUM,   84},
    {Element::ASTATINE,  85}, {Element::RADON,       86}, {Element::FRANCIUM,   87},
    {Element::RADIUM,    88}, {Element::ACTINIUM,    89}, {Element::THORIUM,    90}
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
    {"HOH", Residue::HOH}, {"NA",  Residue::NA},  {"MG",  Residue::MG},  {"CA",  Residue::CA},
    {"ZN",  Residue::ZN},  {"FE",  Residue::FE},  {"FE2", Residue::FE2}, {"CSX", Residue::CSX},
    {"CL",  Residue::CL},

    {"A",   Residue::A},   {"G",   Residue::G},   {"C",   Residue::C},   {"U",   Residue::U},
    {"T",   Residue::T}
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

const std::unordered_map<std::string_view, Remoteness> AtomicInfoHelper::m_remoteness_map
{
    {" ", Remoteness::NONE}, {"A", Remoteness::ALPHA}, {"B", Remoteness::BETA},
    {"G", Remoteness::GAMMA},{"D", Remoteness::DELTA}, {"E", Remoteness::EPSILON},
    {"Z", Remoteness::ZETA}, {"H", Remoteness::ETA},
    {"1", Remoteness::ONE},  {"2", Remoteness::TWO},   {"3", Remoteness::THREE},
    {"4", Remoteness::FOUR}, {"5", Remoteness::FIVE},
    {"X", Remoteness::EXTRA}
};

const std::unordered_map<std::string_view, Branch> AtomicInfoHelper::m_branch_map
{
    {" ", Branch::NONE}, {"1", Branch::ONE}, {"2", Branch::TWO}, {"3", Branch::THREE},
    {"T", Branch::TERMINAL}
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
    {Residue::HOH, "HOH"}, {Residue::NA,  "NA"},  {Residue::MG,  "MG"},  {Residue::CA,  "CA"},
    {Residue::ZN,  "ZN"},  {Residue::FE,  "FE"},  {Residue::FE2, "FE2"}, {Residue::CSX, "CSX"},
    {Residue::CL,  "CL"},
    {Residue::A,   "A"  }, {Residue::G,   "G"  }, {Residue::C,   "C"  }, {Residue::U,   "U"  },
    {Residue::T,   "T"  },
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

const std::unordered_map<Remoteness, std::string> AtomicInfoHelper::m_remoteness_label_map
{
    {Remoteness::NONE, ""}, {Remoteness::ALPHA, "#alpha"},
    {Remoteness::BETA, "#beta"}, {Remoteness::GAMMA, "#gamma"},
    {Remoteness::DELTA, "#delta"}, {Remoteness::EPSILON, "#varepsilon"},
    {Remoteness::ZETA, "#zeta"}, {Remoteness::ETA, "#eta"},
    {Remoteness::ONE, "1"}, {Remoteness::TWO, "2"},
    {Remoteness::THREE, "3"}, {Remoteness::FOUR, "4"},
    {Remoteness::FIVE, "5"}, {Remoteness::EXTRA, "X"},
    {Remoteness::UNK, "UNK"}
};

const std::unordered_map<Remoteness, char> AtomicInfoHelper::m_remoteness_char_map
{
    {Remoteness::NONE, ' '}, {Remoteness::ALPHA, 'A'},
    {Remoteness::BETA, 'B'}, {Remoteness::GAMMA, 'G'},
    {Remoteness::DELTA, 'D'}, {Remoteness::EPSILON, 'E'},
    {Remoteness::ZETA, 'Z'}, {Remoteness::ETA, 'H'},
    {Remoteness::ONE, '1'}, {Remoteness::TWO, '2'},
    {Remoteness::THREE, '3'}, {Remoteness::FOUR, '4'},
    {Remoteness::FIVE, '5'}, {Remoteness::EXTRA, 'X'},
    {Remoteness::UNK, '?'}
};

const std::unordered_map<Branch, std::string> AtomicInfoHelper::m_branch_label_map
{
    {Branch::NONE, ""}, {Branch::ONE, "1"}, {Branch::TWO, "2"}, {Branch::THREE, "3"},
    {Branch::TERMINAL, "T"}, {Branch::UNK, "UNK"}
};

const std::unordered_map<Branch, char> AtomicInfoHelper::m_branch_char_map
{
    {Branch::NONE, ' '}, {Branch::ONE, '1'}, {Branch::TWO, '2'}, {Branch::THREE, '3'},
    {Branch::TERMINAL, 'T'}, {Branch::UNK, '?'}
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

    {Residue::A,     2}, {Residue::G,     3}, {Residue::C,     4}, {Residue::U,     5},
    {Residue::T,     6}
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

    {Residue::A,    2}, {Residue::G,    3}, {Residue::C,    4}, {Residue::U,    5},
    {Residue::T,    6}
};

int AtomicInfoHelper::GetAtomicNumber(Element element)
{
    if (m_atomic_number_map.find(element) == m_atomic_number_map.end()) return 0;
    return m_atomic_number_map.at(element);
}

size_t AtomicInfoHelper::GetGroupClassCount(void)
{
    return m_group_class_key_list.size();
}

size_t AtomicInfoHelper::GetStandardResidueCount(void)
{
    return m_standard_residue_list.size();
}

size_t AtomicInfoHelper::GetElementCount(void)
{
    return m_element_map.size();
}

const std::string & AtomicInfoHelper::GetGroupClassKey(size_t class_id)
{
    if (class_id >= m_group_class_key_list.size())
    {
        throw std::out_of_range("class_id out of range");
    }
    return m_group_class_key_list.at(class_id);
}

const std::string & AtomicInfoHelper::GetElementClassKey(void)
{
    return m_group_class_key_list.at(0);
}

const std::string & AtomicInfoHelper::GetResidueClassKey(void)
{
    return m_group_class_key_list.at(1);
}

const std::string & AtomicInfoHelper::GetStructureClassKey(void)
{
    return m_group_class_key_list.at(2);
}

const std::vector<Residue> & AtomicInfoHelper::GetStandardResidueList(void)
{
    return m_standard_residue_list;
}

const std::vector<Element> & AtomicInfoHelper::GetStandardElementList(void)
{
    return m_standard_element_list;
}

const std::vector<Remoteness> & AtomicInfoHelper::GetStandardRemotenessList(void)
{
    return m_standard_remoteness_list;
}

const std::vector<Branch> & AtomicInfoHelper::GetStandardBranchList(void)
{
    return m_standard_branch_list;
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

Residue AtomicInfoHelper::GetResidueFromString(const std::string & name)
{
    thread_local static std::unordered_map<std::string, int> unknown_name_count_list;
    if (m_residue_map.find(name) == m_residue_map.end())
    {
        if (unknown_name_count_list.find(name) == unknown_name_count_list.end())
        {
            Logger::Log(LogLevel::Warning, 
                        "AtomicInfoHelper::GetResidueFromString - Unknown string: " + name);
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

Remoteness AtomicInfoHelper::GetRemotenessFromString(const std::string & name)
{
    thread_local static std::unordered_map<std::string, int> unknown_name_count_list;
    if (m_remoteness_map.find(name) == m_remoteness_map.end())
    {
        if (unknown_name_count_list.find(name) == unknown_name_count_list.end())
        {
            Logger::Log(LogLevel::Warning, 
                        "AtomicInfoHelper::GetRemotenessFromString - Unknown string: " + name);
            unknown_name_count_list[name] = 1;
        }
        else
        {
            unknown_name_count_list[name]++;
        }
        return Remoteness::UNK;
    }
    return m_remoteness_map.at(name);
}

Branch AtomicInfoHelper::GetBranchFromString(const std::string & name)
{
    thread_local static std::unordered_map<std::string, int> unknown_name_count_list;
    if (m_branch_map.find(name) == m_branch_map.end())
    {
        if (unknown_name_count_list.find(name) == unknown_name_count_list.end())
        {
            Logger::Log(LogLevel::Warning, 
                        "AtomicInfoHelper::GetBranchFromString - Unknown string: " + name);
            unknown_name_count_list[name] = 1;
        }
        else
        {
            unknown_name_count_list[name]++;
        }
        return Branch::UNK;
    }
    return m_branch_map.at(name);
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

const std::string & AtomicInfoHelper::GetLabel(Remoteness remoteness)
{
    static std::string unk_label{"?"};
    if (m_remoteness_label_map.find(remoteness) == m_remoteness_label_map.end())
    {
        Logger::Log(LogLevel::Warning, 
                    "AtomicInfoHelper::GetLabel - Unknown Remoteness: "
                    + std::to_string(static_cast<int>(remoteness)));
        return unk_label;
    }
    return m_remoteness_label_map.at(remoteness);
}

const std::string & AtomicInfoHelper::GetLabel(Branch branch)
{
    static std::string unk_label{"?"};
    if (m_branch_label_map.find(branch) == m_branch_label_map.end())
    {
        Logger::Log(LogLevel::Warning, 
                    "AtomicInfoHelper::GetLabel - Unknown Branch: "
                    + std::to_string(static_cast<int>(branch)));
        return unk_label;
    }
    return m_branch_label_map.at(branch);
}

char AtomicInfoHelper::GetChar(Remoteness remoteness)
{
    static char unk_char{'?'};
    if (m_remoteness_char_map.find(remoteness) == m_remoteness_char_map.end())
    {
        Logger::Log(LogLevel::Warning, 
                    "AtomicInfoHelper::GetChar - Unknown Remoteness: "
                    + std::to_string(static_cast<int>(remoteness)));
        return unk_char;
    }
    return m_remoteness_char_map.at(remoteness);
}

char AtomicInfoHelper::GetChar(Branch branch)
{
    static char unk_char{'?'};
    if (m_branch_char_map.find(branch) == m_branch_char_map.end())
    {
        Logger::Log(LogLevel::Warning, 
                    "AtomicInfoHelper::GetChar - Unknown Branch: "
                    + std::to_string(static_cast<int>(branch)));
        return unk_char;
    }
    return m_branch_char_map.at(branch);
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
