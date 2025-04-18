#include "AtomicInfoHelper.hpp"

const std::vector<int> AtomicInfoHelper::m_standard_residue_list
{
    static_cast<int>(Residue::ALA), static_cast<int>(Residue::ARG),
    static_cast<int>(Residue::ASN), static_cast<int>(Residue::ASP),
    static_cast<int>(Residue::CYS), static_cast<int>(Residue::GLN),
    static_cast<int>(Residue::GLU), static_cast<int>(Residue::GLY),
    static_cast<int>(Residue::HIS), static_cast<int>(Residue::ILE),
    static_cast<int>(Residue::LEU), static_cast<int>(Residue::LYS),
    static_cast<int>(Residue::MET), static_cast<int>(Residue::PHE),
    static_cast<int>(Residue::PRO), static_cast<int>(Residue::SER),
    static_cast<int>(Residue::THR), static_cast<int>(Residue::TRP),
    static_cast<int>(Residue::TYR), static_cast<int>(Residue::VAL)
};

const std::vector<int> AtomicInfoHelper::m_standard_element_list
{
    static_cast<int>(Element::HYDROGEN), static_cast<int>(Element::CARBON),
    static_cast<int>(Element::NITROGEN), static_cast<int>(Element::OXYGEN),
    static_cast<int>(Element::SULFUR)
};

const std::vector<int> AtomicInfoHelper::m_standard_remoteness_list
{
    static_cast<int>(Remoteness::NONE), static_cast<int>(Remoteness::ALPHA),
    static_cast<int>(Remoteness::BETA), static_cast<int>(Remoteness::GAMMA),
    static_cast<int>(Remoteness::DELTA), static_cast<int>(Remoteness::EPSILON),
    static_cast<int>(Remoteness::ZETA), static_cast<int>(Remoteness::ETA)
};

const std::vector<int> AtomicInfoHelper::m_standard_branch_list
{
    static_cast<int>(Branch::NONE), static_cast<int>(Branch::ONE),
    static_cast<int>(Branch::TWO), static_cast<int>(Branch::THREE)
};

const std::unordered_map<Element, int> AtomicInfoHelper::atomic_number_map
{
    {Element::HYDROGEN, 1}, {Element::CARBON,    6}, {Element::NITROGEN,   7},
    {Element::OXYGEN,   8}, {Element::SODIUM,   11}, {Element::MAGNESIUM, 12},
    {Element::SULFUR,  16}, {Element::CHLORINE, 17}, {Element::CALCIUM,   20},
    {Element::IRON,    26}, {Element::ZINC,     30}
};

const std::unordered_map<std::string_view, Residue> AtomicInfoHelper::residue_map
{
    {"ALA", Residue::ALA}, {"ARG", Residue::ARG}, {"ASN", Residue::ASN}, {"ASP", Residue::ASP},
    {"CYS", Residue::CYS}, {"GLN", Residue::GLN}, {"GLU", Residue::GLU}, {"GLY", Residue::GLY},
    {"HIS", Residue::HIS}, {"ILE", Residue::ILE}, {"LEU", Residue::LEU}, {"LYS", Residue::LYS},
    {"MET", Residue::MET}, {"PHE", Residue::PHE}, {"PRO", Residue::PRO}, {"SER", Residue::SER},
    {"THR", Residue::THR}, {"TRP", Residue::TRP}, {"TYR", Residue::TYR}, {"VAL", Residue::VAL},
    {"HOH", Residue::HOH}, {"NA",  Residue::NA},  {"MG",  Residue::MG},  {"CA",  Residue::CA},
    {"ZN",  Residue::ZN},  {"FE",  Residue::FE},  {"FE2", Residue::FE2}, {"CSX", Residue::CSX},
    {"CL",  Residue::CL},  {"UNK", Residue::UNK}
};

const std::unordered_map<std::string_view, Element> AtomicInfoHelper::element_map
{
    {"H",  Element::HYDROGEN},
    {"C",  Element::CARBON},    {"N",  Element::NITROGEN},
    {"O",  Element::OXYGEN},    {"S",  Element::SULFUR},
    {"CA", Element::CALCIUM},   {"ZN", Element::ZINC},    {"NA", Element::SODIUM},
    {"MG", Element::MAGNESIUM}, {"FE", Element::IRON},    {"CL", Element::CHLORINE},
    {"UNK",Element::UNK}
};

const std::unordered_map<std::string_view, Remoteness> AtomicInfoHelper::remoteness_map
{
    {" ", Remoteness::NONE}, {"A", Remoteness::ALPHA}, {"B", Remoteness::BETA},
    {"G", Remoteness::GAMMA},{"D", Remoteness::DELTA}, {"E", Remoteness::EPSILON},
    {"Z", Remoteness::ZETA}, {"H", Remoteness::ETA},
    {"1", Remoteness::ONE},  {"2", Remoteness::TWO},   {"3", Remoteness::THREE},
    {"4", Remoteness::FOUR}, {"5", Remoteness::FIVE},
    {"X", Remoteness::EXTRA},{"UNK",Remoteness::UNK}
};

const std::unordered_map<std::string_view, Branch> AtomicInfoHelper::branch_map
{
    {" ", Branch::NONE}, {"1", Branch::ONE}, {"2", Branch::TWO}, {"3", Branch::THREE},
    {"T", Branch::TERMINAL}, {"UNK", Branch::UNK}
};

const std::unordered_map<Residue, std::string> AtomicInfoHelper::residue_label_map
{
    {Residue::ALA, "ALA"}, {Residue::ARG, "ARG"}, {Residue::ASN, "ASN"}, {Residue::ASP, "ASP"},
    {Residue::CYS, "CYS"}, {Residue::GLN, "GLN"}, {Residue::GLU, "GLU"}, {Residue::GLY, "GLY"},
    {Residue::HIS, "HIS"}, {Residue::ILE, "ILE"}, {Residue::LEU, "LEU"}, {Residue::LYS, "LYS"},
    {Residue::MET, "MET"}, {Residue::PHE, "PHE"}, {Residue::PRO, "PRO"}, {Residue::SER, "SER"},
    {Residue::THR, "THR"}, {Residue::TRP, "TRP"}, {Residue::TYR, "TYR"}, {Residue::VAL, "VAL"},
    {Residue::HOH, "HOH"}, {Residue::NA,  "NA"},  {Residue::MG,  "MG"},  {Residue::CA,  "CA"},
    {Residue::ZN,  "ZN"},  {Residue::FE,  "FE"},  {Residue::FE2, "FE2"}, {Residue::CSX, "CSX"},
    {Residue::CL,  "CL"},  {Residue::UNK, "UNK"}
};

const std::unordered_map<Element, std::string> AtomicInfoHelper::element_label_map
{
    {Element::HYDROGEN, "H"}, {Element::CARBON,    "C"}, {Element::NITROGEN,   "N"},
    {Element::OXYGEN,   "O"}, {Element::SULFUR,    "S"}, {Element::CALCIUM,   "Ca"},
    {Element::ZINC,    "Zn"}, {Element::SODIUM,   "Na"}, {Element::MAGNESIUM, "Mg"},
    {Element::IRON,    "Fe"}, {Element::CHLORINE, "Cl"}, {Element::UNK,        "X"}
};

const std::unordered_map<Remoteness, std::string> AtomicInfoHelper::remoteness_label_map
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

const std::unordered_map<Branch, std::string> AtomicInfoHelper::branch_label_map
{
    {Branch::NONE, ""}, {Branch::ONE, "1"}, {Branch::TWO, "2"}, {Branch::THREE, "3"},
    {Branch::TERMINAL, "T"}, {Branch::UNK, "UNK"}
};

const std::unordered_map<Element, int> AtomicInfoHelper::element_color_map
{
    {Element::HYDROGEN, 921}, {Element::CARBON,   633}, {Element::NITROGEN,  418},
    {Element::OXYGEN,   862}, {Element::SULFUR,   619}, {Element::CALCIUM,     1},
    {Element::ZINC,       1}, {Element::SODIUM,     1}, {Element::MAGNESIUM,   1},
    {Element::IRON,       1}, {Element::CHLORINE,   1}, {Element::UNK,         1}
};