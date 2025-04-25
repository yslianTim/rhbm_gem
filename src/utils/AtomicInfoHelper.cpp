#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"

#include <iostream>
#include <algorithm>

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
    {Element::HYDROGEN, 1}, {Element::CARBON,    6}, {Element::NITROGEN,   7},
    {Element::OXYGEN,   8}, {Element::SODIUM,   11}, {Element::MAGNESIUM, 12},
    {Element::SULFUR,  16}, {Element::CHLORINE, 17}, {Element::CALCIUM,   20},
    {Element::IRON,    26}, {Element::ZINC,     30}
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
    {"CL",  Residue::CL}
};

const std::unordered_map<std::string_view, Element> AtomicInfoHelper::m_element_map
{
    {"H",  Element::HYDROGEN},
    {"C",  Element::CARBON},    {"N",  Element::NITROGEN},
    {"O",  Element::OXYGEN},    {"P",  Element::PHOSPHORUS},
    {"S",  Element::SULFUR},
    {"CA", Element::CALCIUM},   {"ZN", Element::ZINC},    {"NA", Element::SODIUM},
    {"MG", Element::MAGNESIUM}, {"FE", Element::IRON},    {"CL", Element::CHLORINE}
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
    {" ", Structure::FREE}, {"BEND", Structure::BEND}, {"STRN", Structure::STRN},
    {"OTHER", Structure::OTHER},
    {"HELX_P", Structure::HELX_P},
    {"TURN_P", Structure::TURN_P}
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
    {Residue::CL,  "CL"},  {Residue::UNK, "UNK"}
};

const std::unordered_map<Element, std::string> AtomicInfoHelper::m_element_label_map
{
    {Element::HYDROGEN, "H"}, {Element::CARBON,     "C"}, {Element::NITROGEN,   "N"},
    {Element::OXYGEN,   "O"}, {Element::PHOSPHORUS, "P"},
    {Element::SULFUR,   "S"}, {Element::CALCIUM,   "Ca"},
    {Element::ZINC,    "Zn"}, {Element::SODIUM,    "Na"}, {Element::MAGNESIUM, "Mg"},
    {Element::IRON,    "Fe"}, {Element::CHLORINE,  "Cl"}, {Element::UNK,        "X"}
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

const std::unordered_map<Branch, std::string> AtomicInfoHelper::m_branch_label_map
{
    {Branch::NONE, ""}, {Branch::ONE, "1"}, {Branch::TWO, "2"}, {Branch::THREE, "3"},
    {Branch::TERMINAL, "T"}, {Branch::UNK, "UNK"}
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
    {Residue::THR, 880}, {Residue::TRP, 881}, {Residue::TYR, 882}, {Residue::VAL, 883}
};

const std::unordered_map<Element, int> AtomicInfoHelper::m_element_marker_map
{
    {Element::HYDROGEN,  5}, {Element::CARBON,    53}, {Element::NITROGEN,  55},
    {Element::OXYGEN,   59}, {Element::PHOSPHORUS,44},
    {Element::SULFUR,   27}, {Element::CALCIUM,   28},
    {Element::ZINC,     30}, {Element::SODIUM,    31}, {Element::MAGNESIUM, 32},
    {Element::IRON,     35}, {Element::CHLORINE,  36}, {Element::UNK,       37}
};

const std::unordered_map<Residue, int> AtomicInfoHelper::m_residue_marker_map
{
    {Residue::ALA, 20}, {Residue::ARG, 21}, {Residue::ASN, 22}, {Residue::ASP, 23},
    {Residue::CYS, 24}, {Residue::GLN, 25}, {Residue::GLU, 26}, {Residue::GLY, 27},
    {Residue::HIS, 28}, {Residue::ILE, 29}, {Residue::LEU, 30}, {Residue::LYS, 31},
    {Residue::MET, 32}, {Residue::PHE, 33}, {Residue::PRO, 34}, {Residue::SER, 35},
    {Residue::THR, 36}, {Residue::TRP, 37}, {Residue::TYR, 38}, {Residue::VAL, 39}
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

bool AtomicInfoHelper::IsStandardResidue(Residue residue)
{
    return std::find(m_standard_residue_list.begin(), m_standard_residue_list.end(), residue)
            != m_standard_residue_list.end();
}

Residue AtomicInfoHelper::GetResidueFromString(const std::string & name)
{
    if (m_residue_map.find(name) == m_residue_map.end()) return Residue::UNK;
    return m_residue_map.at(name);
}

Element AtomicInfoHelper::GetElementFromString(const std::string & name)
{
    if (m_element_map.find(name) == m_element_map.end()) return Element::UNK;
    return m_element_map.at(name);
}

Remoteness AtomicInfoHelper::GetRemotenessFromString(const std::string & name)
{
    if (m_remoteness_map.find(name) == m_remoteness_map.end()) return Remoteness::UNK;
    return m_remoteness_map.at(name);
}

Branch AtomicInfoHelper::GetBranchFromString(const std::string & name)
{
    if (m_branch_map.find(name) == m_branch_map.end()) return Branch::UNK;
    return m_branch_map.at(name);
}

Structure AtomicInfoHelper::GetStructureFromString(const std::string & name)
{
    if (m_structure_map.find(name) == m_structure_map.end()) return Structure::UNK;
    return m_structure_map.at(name);
}

const std::string & AtomicInfoHelper::GetLabel(Residue residue)
{
    return m_residue_label_map.at(residue);
}

const std::string & AtomicInfoHelper::GetLabel(Element element)
{
    return m_element_label_map.at(element);
}

const std::string & AtomicInfoHelper::GetLabel(Remoteness remoteness)
{
    return m_remoteness_label_map.at(remoteness);
}

const std::string & AtomicInfoHelper::GetLabel(Branch branch)
{
    return m_branch_label_map.at(branch);
}

short AtomicInfoHelper::GetDisplayColor(Element element)
{
    return static_cast<short>(m_element_color_map.at(element));
}

short AtomicInfoHelper::GetDisplayColor(Residue residue)
{
    return static_cast<short>(m_residue_color_map.at(residue));
}

short AtomicInfoHelper::GetDisplayMarker(Element element)
{
    return static_cast<short>(m_element_marker_map.at(element));
}

short AtomicInfoHelper::GetDisplayMarker(Residue residue)
{
    return static_cast<short>(m_residue_marker_map.at(residue));
}