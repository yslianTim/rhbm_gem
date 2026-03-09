#include <rhbm_gem/utils/ComponentHelper.hpp>
#include <rhbm_gem/utils/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/Logger.hpp>

#include <stdexcept>
#include <mutex>

const std::unordered_map<Residue, std::vector<Spot>> ComponentHelper::m_spot_map
{
    {Residue::ALA, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB}},
    {Residue::ARG, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::CD,  Spot::NE,  Spot::CZ,  Spot::NH1,
                    Spot::NH2}},
    {Residue::ASN, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::OD1, Spot::ND2}},
    {Residue::ASP, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::OD1, Spot::OD2}},
    {Residue::CYS, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::SG}},
    {Residue::GLN, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::CD,  Spot::OE1, Spot::NE2}},
    {Residue::GLU, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::CD,  Spot::OE1, Spot::OE2}},
    {Residue::GLY, {Spot::C,   Spot::CA,  Spot::N,   Spot::O}},
    {Residue::HIS, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::ND1, Spot::CD2, Spot::CE1, Spot::NE2}},
    {Residue::ILE, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG1, Spot::CG2, Spot::CD1}},
    {Residue::LEU, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::CD1, Spot::CD2}},
    {Residue::LYS, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::CD,  Spot::CE,  Spot::NZ}},
    {Residue::MET, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::SD,  Spot::CE}},
    {Residue::PHE, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::CD1, Spot::CD2, Spot::CE1, Spot::CE2,
                    Spot::CZ}},
    {Residue::PRO, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::CD}},
    {Residue::SER, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::OG}},
    {Residue::THR, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::OG1, Spot::CG2}},
    {Residue::TRP, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::CD1, Spot::CD2, Spot::NE1, Spot::CE2,
                    Spot::CE3, Spot::CZ2, Spot::CZ3, Spot::CH2}},
    {Residue::TYR, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG,  Spot::CD1, Spot::CD2, Spot::CE1, Spot::CE2,
                    Spot::CZ,  Spot::OH}},
    {Residue::VAL, {Spot::C,   Spot::CA,  Spot::N,   Spot::O,   Spot::CB,
                    Spot::CG1, Spot::CG2}},
    
    {Residue::A,   {Spot::P,   Spot::OP1, Spot::OP2,
                    Spot::O5p, Spot::C5p, Spot::C4p, Spot::O4p, Spot::C3p,
                    Spot::O3p, Spot::C2p, Spot::O2p, Spot::C1p,
                    Spot::N9,  Spot::C8,  Spot::N7,  Spot::C5,  Spot::C6,
                    Spot::N6,  Spot::N1,  Spot::C2,  Spot::N3,  Spot::C4}},
    {Residue::C,   {Spot::P,   Spot::OP1, Spot::OP2,
                    Spot::O5p, Spot::C5p, Spot::C4p, Spot::O4p, Spot::C3p,
                    Spot::O3p, Spot::C2p, Spot::O2p, Spot::C1p,
                    Spot::N1,  Spot::C2,  Spot::O2,  Spot::N3,  Spot::C4,
                    Spot::N4,  Spot::C5,  Spot::C6}},
    {Residue::G,   {Spot::P,   Spot::OP1, Spot::OP2,
                    Spot::O5p, Spot::C5p, Spot::C4p, Spot::O4p, Spot::C3p,
                    Spot::O3p, Spot::C2p, Spot::O2p, Spot::C1p,
                    Spot::N9,  Spot::C8,  Spot::N7,  Spot::C5,  Spot::C6,
                    Spot::O6,  Spot::N1,  Spot::C2,  Spot::N2,  Spot::N3,
                    Spot::C4}},
    {Residue::U,   {Spot::P,   Spot::OP1, Spot::OP2,
                    Spot::O5p, Spot::C5p, Spot::C4p, Spot::O4p, Spot::C3p,
                    Spot::O3p, Spot::C2p, Spot::O2p, Spot::C1p,
                    Spot::N1,  Spot::C2,  Spot::O2,  Spot::N3,  Spot::C4,
                    Spot::O4,  Spot::C5,  Spot::C6}},
    {Residue::DA,  {Spot::P,   Spot::OP1, Spot::OP2,
                    Spot::O5p, Spot::C5p, Spot::C4p, Spot::O4p, Spot::C3p,
                    Spot::O3p, Spot::C2p, Spot::C1p,
                    Spot::N9,  Spot::C8,  Spot::N7,  Spot::C5,  Spot::C6,
                    Spot::N6,  Spot::N1,  Spot::C2,  Spot::N3,  Spot::C4}},
    {Residue::DC,  {Spot::P,   Spot::OP1, Spot::OP2,
                    Spot::O5p, Spot::C5p, Spot::C4p, Spot::O4p, Spot::C3p,
                    Spot::O3p, Spot::C2p, Spot::C1p,
                    Spot::N1,  Spot::C2,  Spot::O2,  Spot::N3,  Spot::C4,
                    Spot::N4,  Spot::C5,  Spot::C6}},
    {Residue::DG,  {Spot::P,   Spot::OP1, Spot::OP2,
                    Spot::O5p, Spot::C5p, Spot::C4p, Spot::O4p, Spot::C3p,
                    Spot::O3p, Spot::C2p, Spot::C1p,
                    Spot::N9,  Spot::C8,  Spot::N7,  Spot::C5,  Spot::C6,
                    Spot::O6,  Spot::N1,  Spot::C2,  Spot::N2,  Spot::N3,
                    Spot::C4}},
    {Residue::DT,  {Spot::P,   Spot::OP1, Spot::OP2,
                    Spot::O5p, Spot::C5p, Spot::C4p, Spot::O4p, Spot::C3p,
                    Spot::O3p, Spot::C2p, Spot::C1p,
                    Spot::N1,  Spot::C2,  Spot::O2,  Spot::N3,  Spot::C4,
                    Spot::O4,  Spot::C5,  Spot::C7,  Spot::C6}}
};

const std::unordered_map<Residue, std::vector<Link>> ComponentHelper::m_link_map
{
    {Residue::ALA, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB}},
    {Residue::ARG, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_CD,   Link::CD_NE,   Link::NE_CZ,   Link::CZ_NH1,
                    Link::CZ_NH2}},
    {Residue::ASN, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_OD1,  Link::CG_ND2}},
    {Residue::ASP, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_OD1,  Link::CG_OD2}},
    {Residue::CYS, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_SG}},
    {Residue::GLN, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_CD,   Link::CD_OE1,  Link::CD_NE2}},
    {Residue::GLU, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_CD,   Link::CD_OE1,  Link::CD_OE2}},
    {Residue::GLY, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N}},
    {Residue::HIS, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_ND1,  Link::CG_CD2,  Link::ND1_CE1, Link::CD2_NE2,
                    Link::CE1_NE2}},
    {Residue::ILE, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG1,  Link::CB_CG2,  Link::CG1_CD1}},
    {Residue::LEU, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_CD1,  Link::CG_CD2}},
    {Residue::LYS, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_CD,   Link::CD_CE,   Link::CE_NZ}},
    {Residue::MET, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_SD,   Link::SD_CE}},
    {Residue::PHE, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_CD1,  Link::CG_CD2,  Link::CD1_CE1, Link::CD2_CE2,
                    Link::CE1_CZ,  Link::CE2_CZ}},
    {Residue::PRO, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_CD,   Link::N_CD}},
    {Residue::SER, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_OG}},
    {Residue::THR, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_OG1,  Link::CB_CG2}},
    {Residue::TRP, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_CD1,  Link::CG_CD2,  Link::CD1_NE1, Link::CD2_CE2,
                    Link::CD2_CE3, Link::NE1_CE2, Link::CE2_CZ2, Link::CE3_CZ3, Link::CZ2_CH2,
                    Link::CZ3_CH2}},
    {Residue::TYR, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG,   Link::CG_CD1,  Link::CG_CD2,  Link::CD1_CE1, Link::CD2_CE2,
                    Link::CE1_CZ,  Link::CE2_CZ,  Link::CZ_OH}},
    {Residue::VAL, {Link::N_CA,    Link::CA_C,    Link::C_O,     Link::C_N,     Link::CA_CB,
                    Link::CB_CG1,  Link::CB_CG2}},
    
    {Residue::A,   {Link::P_OP1,   Link::P_OP2,   Link::P_O5p,
                    Link::O5p_C5p, Link::C5p_C4p, Link::C4p_O4p, Link::C4p_C3p,
                    Link::O4p_C1p, Link::C3p_O3p, Link::C3p_C2p, Link::C2p_O2p, Link::C2p_C1p,
                    Link::C1p_N9,  Link::N9_C8,   Link::N9_C4,   Link::C8_N7,   Link::N7_C5,
                    Link::C5_C6,   Link::C5_C4,   Link::C6_N6,   Link::C6_N1,   Link::N1_C2,
                    Link::C2_N3,   Link::N3_C4}},
    {Residue::C,   {Link::P_OP1,   Link::P_OP2,   Link::P_O5p,
                    Link::O5p_C5p, Link::C5p_C4p, Link::C4p_O4p, Link::C4p_C3p,
                    Link::O4p_C1p, Link::C3p_O3p, Link::C3p_C2p, Link::C2p_O2p, Link::C2p_C1p,
                    Link::C1p_N1,  Link::N1_C2,   Link::N1_C6,   Link::C2_O2,   Link::C2_N3,
                    Link::N3_C4,   Link::C4_N4,   Link::C4_C5,   Link::C5_C6}},
    {Residue::G,   {Link::P_OP1,   Link::P_OP2,   Link::P_O5p,
                    Link::O5p_C5p, Link::C5p_C4p, Link::C4p_O4p, Link::C4p_C3p,
                    Link::O4p_C1p, Link::C3p_O3p, Link::C3p_C2p, Link::C2p_O2p, Link::C2p_C1p,
                    Link::C1p_N9,  Link::N9_C8,   Link::N9_C4,   Link::C8_N7,   Link::N7_C5,
                    Link::C5_C6,   Link::C5_C4,   Link::C6_O6,   Link::C6_N1,   Link::N1_C2,
                    Link::C2_N2,   Link::C2_N3,   Link::N3_C4}},
    {Residue::U,   {Link::P_OP1,   Link::P_OP2,   Link::P_O5p,
                    Link::O5p_C5p, Link::C5p_C4p, Link::C4p_O4p, Link::C4p_C3p,
                    Link::O4p_C1p, Link::C3p_O3p, Link::C3p_C2p, Link::C2p_O2p, Link::C2p_C1p,
                    Link::C1p_N1,  Link::N1_C2,   Link::N1_C6,   Link::C2_O2,   Link::C2_N3,
                    Link::N3_C4,   Link::C4_O4,   Link::C4_C5,   Link::C5_C6}},
    {Residue::DA,  {Link::P_OP1,   Link::P_OP2,   Link::P_O5p,
                    Link::O5p_C5p, Link::C5p_C4p, Link::C4p_O4p, Link::C4p_C3p,
                    Link::O4p_C1p, Link::C3p_O3p, Link::C3p_C2p, Link::C2p_C1p,
                    Link::C1p_N9,  Link::N9_C8,   Link::N9_C4,   Link::C8_N7,   Link::N7_C5,
                    Link::C5_C6,   Link::C5_C4,   Link::C6_N6,   Link::C6_N1,   Link::N1_C2,
                    Link::C2_N3,   Link::N3_C4}},
    {Residue::DC,  {Link::P_OP1,   Link::P_OP2,   Link::P_O5p,
                    Link::O5p_C5p, Link::C5p_C4p, Link::C4p_O4p, Link::C4p_C3p,
                    Link::O4p_C1p, Link::C3p_O3p, Link::C3p_C2p, Link::C2p_C1p,
                    Link::C1p_N1,  Link::N1_C2,   Link::N1_C6,   Link::C2_O2,   Link::C2_N3,
                    Link::N3_C4,   Link::C4_N4,   Link::C4_C5,   Link::C5_C6}},
    {Residue::DG,  {Link::P_OP1,   Link::P_OP2,   Link::P_O5p,
                    Link::O5p_C5p, Link::C5p_C4p, Link::C4p_O4p, Link::C4p_C3p,
                    Link::O4p_C1p, Link::C3p_O3p, Link::C3p_C2p, Link::C2p_C1p,
                    Link::C1p_N9,  Link::N9_C8,   Link::N9_C4,   Link::C8_N7,   Link::N7_C5,
                    Link::C5_C6,   Link::C5_C4,   Link::C6_O6,   Link::C6_N1,   Link::N1_C2,
                    Link::C2_N2,   Link::C2_N3,   Link::N3_C4}},
    {Residue::DT,  {Link::P_OP1,   Link::P_OP2,   Link::P_O5p,
                    Link::O5p_C5p, Link::C5p_C4p, Link::C4p_O4p, Link::C4p_C3p,
                    Link::O4p_C1p, Link::C3p_O3p, Link::C3p_C2p, Link::C2p_C1p,
                    Link::C1p_N1,  Link::N1_C2,   Link::N1_C6,   Link::C2_O2,   Link::C2_N3,
                    Link::N3_C4,   Link::C4_O4,   Link::C4_C5,   Link::C5_C7,   Link::C5_C6}}
};

const std::unordered_map<Residue, std::vector<Element>> ComponentHelper::m_element_map
{
    {Residue::ALA, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON}},

    {Residue::ARG, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON, Element::NITROGEN,
                    Element::CARBON, Element::NITROGEN, Element::NITROGEN}},

    {Residue::ASN, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::OXYGEN, Element::NITROGEN}},

    {Residue::ASP, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::OXYGEN, Element::OXYGEN}},

    {Residue::CYS, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::SULFUR}},

    {Residue::GLN, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON, Element::OXYGEN,
                    Element::NITROGEN}},

    {Residue::GLU, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON, Element::OXYGEN,
                    Element::OXYGEN}},

    {Residue::GLY, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN}},

    {Residue::HIS, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::NITROGEN, Element::CARBON,
                    Element::CARBON, Element::NITROGEN}},

    {Residue::ILE, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON, Element::CARBON}},

    {Residue::LEU, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON, Element::CARBON}},

    {Residue::LYS, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON, Element::CARBON,
                    Element::NITROGEN}},

    {Residue::MET, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::SULFUR, Element::CARBON}},

    {Residue::PHE, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON, Element::CARBON,
                    Element::CARBON, Element::CARBON, Element::CARBON}},

    {Residue::PRO, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON}},

    {Residue::SER, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::OXYGEN}},

    {Residue::THR, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::OXYGEN, Element::CARBON}},

    {Residue::TRP, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::CARBON, Element::CARBON,
                    Element::CARBON, Element::CARBON}},

    {Residue::TYR, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON, Element::CARBON,
                    Element::CARBON, Element::CARBON, Element::CARBON, Element::OXYGEN}},

    {Residue::VAL, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::CARBON, Element::CARBON}},
    
    {Residue::A,   {Element::PHOSPHORUS, Element::OXYGEN,   Element::OXYGEN,
                    Element::OXYGEN, Element::CARBON, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::OXYGEN, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::NITROGEN, Element::CARBON, Element::CARBON,
                    Element::NITROGEN, Element::NITROGEN, Element::CARBON, Element::NITROGEN, Element::CARBON}},
    {Residue::C,   {Element::PHOSPHORUS, Element::OXYGEN,   Element::OXYGEN,
                    Element::OXYGEN, Element::CARBON, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::OXYGEN, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::OXYGEN, Element::NITROGEN, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::CARBON}},
    {Residue::G,   {Element::PHOSPHORUS, Element::OXYGEN,   Element::OXYGEN,
                    Element::OXYGEN, Element::CARBON, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::OXYGEN, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::NITROGEN, Element::CARBON, Element::CARBON,
                    Element::OXYGEN, Element::NITROGEN, Element::CARBON, Element::NITROGEN, Element::NITROGEN,
                    Element::CARBON}},
    {Residue::U,   {Element::PHOSPHORUS, Element::OXYGEN,   Element::OXYGEN,
                    Element::OXYGEN, Element::CARBON, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::OXYGEN, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::OXYGEN, Element::NITROGEN, Element::CARBON,
                    Element::OXYGEN, Element::CARBON, Element::CARBON}},
    {Residue::DA,  {Element::PHOSPHORUS, Element::OXYGEN,   Element::OXYGEN,
                    Element::OXYGEN, Element::CARBON, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::OXYGEN, Element::CARBON, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::NITROGEN, Element::CARBON, Element::CARBON,
                    Element::NITROGEN, Element::NITROGEN, Element::CARBON, Element::NITROGEN, Element::CARBON}},
    {Residue::DC,  {Element::PHOSPHORUS, Element::OXYGEN,   Element::OXYGEN,
                    Element::OXYGEN, Element::CARBON, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::OXYGEN, Element::CARBON, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::OXYGEN, Element::NITROGEN, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::CARBON}},
    {Residue::DG,  {Element::PHOSPHORUS, Element::OXYGEN,   Element::OXYGEN,
                    Element::OXYGEN, Element::CARBON, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::OXYGEN, Element::CARBON, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::NITROGEN, Element::CARBON, Element::CARBON,
                    Element::OXYGEN, Element::NITROGEN, Element::CARBON, Element::NITROGEN, Element::NITROGEN,
                    Element::CARBON}},
    {Residue::DT,  {Element::PHOSPHORUS, Element::OXYGEN,   Element::OXYGEN,
                    Element::OXYGEN, Element::CARBON, Element::CARBON, Element::OXYGEN, Element::CARBON,
                    Element::OXYGEN, Element::CARBON, Element::CARBON,
                    Element::NITROGEN, Element::CARBON, Element::OXYGEN, Element::NITROGEN, Element::CARBON,
                    Element::OXYGEN, Element::CARBON, Element::CARBON, Element::CARBON}}
};

const std::unordered_map<Residue, std::vector<double>> ComponentHelper::m_amber95_partial_charge_map
{
    //               C      CA     N      O      CB
    {Residue::ALA, { 0.597, 0.034,-0.416,-0.568,-0.183}},
    {Residue::ARG, { 0.734,-0.264,-0.348,-0.589,-0.001, 0.039, 0.049,-0.530, 0.808,-0.823,-0.823}},
    {Residue::ASN, { 0.597, 0.014,-0.416,-0.568,-0.204, 0.713,-0.593,-0.919}},
    {Residue::ASP, { 0.537, 0.038,-0.516,-0.582,-0.030, 0.799,-0.801,-0.801}},
    {Residue::CYS, { 0.597, 0.021,-0.416,-0.568,-0.123,-0.285}},
    {Residue::GLN, { 0.597,-0.003,-0.416,-0.568,-0.004,-0.065, 0.695,-0.609,-0.941}},
    {Residue::GLU, { 0.537, 0.040,-0.416,-0.582, 0.056, 0.014, 0.805,-0.819,-0.819}},
    {Residue::GLY, { 0.597,-0.025,-0.416,-0.568}},
    {Residue::HIS, { 0.597,-0.058,-0.416,-0.568,-0.007, 0.187,-0.543,-0.221, 0.164,-0.280}},
    {Residue::ILE, { 0.597, 0.060,-0.416,-0.568, 0.130,-0.043,-0.320,-0.066}},
    {Residue::LEU, { 0.597,-0.052,-0.416,-0.568,-0.110, 0.353,-0.412,-0.412}},
    {Residue::LYS, { 0.734,-0.240,-0.416,-0.589,-0.009, 0.019,-0.048,-0.014,-0.385}},
    {Residue::MET, { 0.597,-0.024,-0.416,-0.568, 0.003, 0.002,-0.054,-0.274}},
    {Residue::PHE, { 0.597,-0.002,-0.416,-0.568,-0.034, 0.012,-0.126,-0.126,-0.170,-0.170,-0.107}},
    {Residue::PRO, { 0.590,-0.027,-0.255,-0.548,-0.007, 0.019, 0.019}},
    {Residue::SER, { 0.597,-0.025,-0.416,-0.568, 0.212,-0.655}},
    {Residue::THR, { 0.597,-0.039,-0.416,-0.568, 0.365,-0.676,-0.244}},
    {Residue::TRP, { 0.597,-0.028,-0.416,-0.568,-0.005,-0.142,-0.164, 0.124,-0.342, 0.138,-0.239,-0.260,-0.197,-0.113}},
    {Residue::TYR, { 0.597,-0.001,-0.416,-0.568,-0.015,-0.001,-0.191,-0.191,-0.234,-0.234, 0.326,-0.558}},
    {Residue::VAL, { 0.597,-0.088,-0.416,-0.568, 0.299,-0.319,-0.319}}
};

const std::unordered_map<Residue, std::vector<double>> ComponentHelper::m_buried_partial_charge_map
{
    //               C      CA     N      O      CB
    {Residue::ALA, { 0.560, 0.267,-0.592,-0.673,-0.392}},
    {Residue::ARG, { 0.558, 0.240,-0.590,-0.658,-0.215,-0.190, 0.116,-0.604, 0.901,-0.908,-0.906}},
    {Residue::ASN, { 0.561, 0.244,-0.585,-0.664,-0.238, 0.673,-0.659,-0.872}},
    {Residue::ASP, { 0.570, 0.246,-0.568,-0.673,-0.276, 0.704,-0.771,-0.761}},
    {Residue::CYS, { 0.574, 0.186,-0.558,-0.694, 0.096,-0.285}},
    {Residue::GLN, { 0.562, 0.234,-0.584,-0.666,-0.212,-0.191, 0.667,-0.653,-0.877}},
    {Residue::GLU, { 0.568, 0.234,-0.575,-0.681,-0.207,-0.223, 0.693,-0.772,-0.764}},
    {Residue::GLY, { 0.581, 0.108,-0.566,-0.668}},
    {Residue::HIS, { 0.563, 0.235,-0.584,-0.663,-0.153, 0.022,-0.296,-0.182,-0.012,-0.225}},
    {Residue::ILE, { 0.560, 0.214,-0.590,-0.657,-0.032,-0.154,-0.360,-0.332}},
    {Residue::LEU, { 0.557, 0.242,-0.593,-0.663,-0.235, 0.049,-0.359,-0.361}},
    {Residue::LYS, { 0.561, 0.238,-0.585,-0.652,-0.214,-0.164,-0.191, 0.037,-0.398}},
    {Residue::MET, { 0.559, 0.244,-0.592,-0.665,-0.274, 0.073,-0.141,-0.070}},
    {Residue::PHE, { 0.559, 0.235,-0.589,-0.656,-0.156, 0.007,-0.204,-0.205,-0.177,-0.179,-0.182}},
    {Residue::PRO, { 0.567, 0.170,-0.310,-0.670,-0.196,-0.165, 0.021}},
    {Residue::SER, { 0.567, 0.193,-0.579,-0.671, 0.079,-0.489}},
    {Residue::THR, { 0.565, 0.169,-0.573,-0.659, 0.232,-0.494,-0.413}},
    {Residue::TRP, { 0.565, 0.233,-0.589,-0.669,-0.133,-0.175,-0.083, 0.013,-0.310, 0.195,-0.230,-0.268,-0.213,-0.190}},
    {Residue::TYR, { 0.561, 0.235,-0.586,-0.662,-0.148,-0.044,-0.182,-0.182,-0.334,-0.345, 0.403,-0.568}},
    {Residue::VAL, { 0.560, 0.215,-0.588,-0.661,-0.005,-0.362,-0.360}}
};

const std::unordered_map<Residue, std::vector<double>> ComponentHelper::m_helix_partial_charge_map
{
    //               C      CA     N      O      CB
    {Residue::ALA, { 0.559, 0.272,-0.598,-0.701,-0.400}},
    {Residue::ARG, { 0.560, 0.243,-0.595,-0.682,-0.224,-0.193, 0.113,-0.596, 0.902,-0.909,-0.911}},
    {Residue::ASN, { 0.559, 0.243,-0.589,-0.693,-0.247, 0.680,-0.659,-0.877}},
    {Residue::ASP, { 0.573, 0.246,-0.574,-0.698,-0.285, 0.713,-0.777,-0.782}},
    {Residue::CYS, { 0.571, 0.188,-0.577,-0.718, 0.086,-0.285}},
    {Residue::GLN, { 0.565, 0.236,-0.587,-0.687,-0.221,-0.193, 0.667,-0.646,-0.885}},
    {Residue::GLU, { 0.569, 0.236,-0.584,-0.705,-0.214,-0.224, 0.694,-0.777,-0.771}},
    {Residue::GLY, { 0.575, 0.119,-0.575,-0.706}},
    {Residue::HIS, { 0.571, 0.235,-0.589,-0.690,-0.167, 0.032,-0.280,-0.185,-0.015,-0.236}},
    {Residue::ILE, { 0.562, 0.217,-0.601,-0.698,-0.039,-0.153,-0.360,-0.331}},
    {Residue::LEU, { 0.558, 0.245,-0.600,-0.699,-0.245, 0.047,-0.359,-0.362}},
    {Residue::LYS, { 0.560, 0.242,-0.592,-0.683,-0.224,-0.166,-0.192, 0.038,-0.402}},
    {Residue::MET, { 0.559, 0.249,-0.602,-0.706,-0.283, 0.072,-0.155,-0.071}},
    {Residue::PHE, { 0.569, 0.239,-0.603,-0.688,-0.171, 0.011,-0.204,-0.204,-0.181,-0.181,-0.185}},
    {Residue::PRO, { 0.555, 0.166,-0.343,-0.727,-0.207,-0.165, 0.021}}, // same as buried charge
    {Residue::SER, { 0.562, 0.197,-0.595,-0.698, 0.071,-0.501}},
    {Residue::THR, { 0.563, 0.173,-0.589,-0.690, 0.225,-0.506,-0.410}},
    {Residue::TRP, { 0.568, 0.234,-0.602,-0.697,-0.142,-0.171,-0.085, 0.018,-0.310, 0.197,-0.225,-0.269,-0.213,-0.188}},
    {Residue::TYR, { 0.573, 0.238,-0.599,-0.690,-0.165,-0.037,-0.181,-0.182,-0.335,-0.349, 0.401,-0.565}},
    {Residue::VAL, { 0.559, 0.220,-0.602,-0.695,-0.013,-0.362,-0.362}}
};

const std::unordered_map<Residue, std::vector<double>> ComponentHelper::m_sheet_partial_charge_map
{
    //               C      CA     N      O      CB
    {Residue::ALA, { 0.562, 0.264,-0.588,-0.661,-0.380}},
    {Residue::ARG, { 0.557, 0.238,-0.586,-0.658,-0.206,-0.189, 0.119,-0.610, 0.901,-0.905,-0.898}},
    {Residue::ASN, { 0.564, 0.239,-0.576,-0.656,-0.233, 0.672,-0.643,-0.871}},
    {Residue::ASP, { 0.570, 0.242,-0.560,-0.665,-0.271, 0.701,-0.769,-0.749}},
    {Residue::CYS, { 0.584, 0.183,-0.550,-0.695, 0.103,-0.285}},
    {Residue::GLN, { 0.563, 0.233,-0.578,-0.673,-0.203,-0.189, 0.664,-0.652,-0.869}},
    {Residue::GLU, { 0.569, 0.228,-0.564,-0.676,-0.199,-0.225, 0.691,-0.763,-0.757}},
    {Residue::GLY, { 0.580, 0.103,-0.563,-0.653}},
    {Residue::HIS, { 0.564, 0.234,-0.577,-0.657,-0.145, 0.020,-0.295,-0.180,-0.011,-0.228}},
    {Residue::ILE, { 0.565, 0.213,-0.587,-0.655,-0.027,-0.153,-0.359,-0.331}},
    {Residue::LEU, { 0.563, 0.239,-0.586,-0.649,-0.227, 0.049,-0.358,-0.360}},
    {Residue::LYS, { 0.562, 0.237,-0.586,-0.662,-0.204,-0.162,-0.191, 0.040,-0.398}},
    {Residue::MET, { 0.563, 0.240,-0.583,-0.654,-0.266, 0.074,-0.134,-0.068}},
    {Residue::PHE, { 0.559, 0.234,-0.581,-0.644,-0.148, 0.008,-0.200,-0.205,-0.173,-0.177,-0.180}},
    {Residue::PRO, { 0.555, 0.166,-0.343,-0.727,-0.207,-0.165, 0.021}}, // same as buried charge
    {Residue::SER, { 0.566, 0.190,-0.571,-0.668, 0.085,-0.477}},
    {Residue::THR, { 0.568, 0.168,-0.569,-0.663, 0.238,-0.482,-0.409}},
    {Residue::TRP, { 0.562, 0.229,-0.575,-0.661,-0.124,-0.181,-0.091, 0.008,-0.311, 0.195,-0.234,-0.270,-0.206,-0.186}},
    {Residue::TYR, { 0.561, 0.233,-0.578,-0.663,-0.138,-0.047,-0.180,-0.183,-0.330,-0.344, 0.408,-0.560}},
    {Residue::VAL, { 0.563, 0.214,-0.587,-0.657, 0.000,-0.361,-0.356}}
};

size_t ComponentHelper::GetAtomCount(Residue residue)
{
    if (m_spot_map.find(residue) == m_spot_map.end())
    {
        throw std::out_of_range(
            "ComponentHelper::GetAtomCount - residue is not supported");
    }
    return m_spot_map.at(residue).size();
}

size_t ComponentHelper::GetBondCount(Residue residue)
{
    if (m_link_map.find(residue) == m_link_map.end())
    {
        Logger::Log(LogLevel::Warning,
            "ComponentHelper::GetBondCount ‑ residue is not supported");
        return 0;
    }
    return m_link_map.at(residue).size();
}

double ComponentHelper::GetPartialCharge(
    Residue residue, Spot spot, Structure structure,
    bool use_amber_table, bool verbose)
{
    // one cache bucket per residue, lazily initialised on first use
    static std::unordered_map<Residue, std::unordered_map<Spot, double>> cache;
    static std::mutex cache_mutex;

    std::lock_guard<std::mutex> lock(cache_mutex);
    auto & residue_cache{ cache[residue] };
    try
    {
        if (residue_cache.empty()) // first request for this residue
        {
            auto atom_size{ m_spot_map.at(residue).size() };
            const auto & spot_list{ m_spot_map.at(residue) };
            const auto & charge_list
            {
                (use_amber_table == true) ?
                GetPartialChargeListAmber(residue) :
                GetPartialChargeList(residue, structure)
            };

            // the four vectors should guaranteed aligned
            if (atom_size != charge_list.size())
            {
                throw std::range_error(
                    "ComponentHelper::GetPartialCharge ‑ the four vectors are not aligned");
            }

            residue_cache.reserve(atom_size);
            for (std::size_t i = 0; i < atom_size; i++)
            {
                residue_cache.emplace(spot_list[i], charge_list[i]);
            }
        }
    }
    catch(const std::exception & except)
    {
        if (verbose == true)
        {
            Logger::Log(LogLevel::Warning,
                "ComponentHelper::GetPartialCharge ‑ " + std::string(except.what()));
        }
        return 0.0;
    }
    
    if (residue_cache.find(spot) != residue_cache.end())
    {
        return residue_cache.at(spot);
    }

    if (verbose == true)
    {
        Logger::Log(LogLevel::Warning, "No partial charge data for this atom.");
    }
    return 0.0;
}

const std::vector<double> & ComponentHelper::GetPartialChargeList(
    Residue residue, Structure structure)
{
    if (structure == Structure::FREE)
    {
        return m_buried_partial_charge_map.at(residue);
    }
    else if (structure == Structure::SHEET)
    {
        return m_sheet_partial_charge_map.at(residue);
    }
    else if (structure >= Structure::HELX_P && structure < Structure::TURN_P)
    {
        return m_helix_partial_charge_map.at(residue);
    }
    else
    {
        throw std::out_of_range(
            "ComponentHelper::GetPartialChargeList ‑ structure is not supported");
    }
}

const std::vector<double> & ComponentHelper::GetPartialChargeListAmber(Residue residue)
{
    return m_amber95_partial_charge_map.at(residue);
}

const std::vector<Spot> & ComponentHelper::GetSpotList(Residue residue)
{
    if (m_spot_map.find(residue) == m_spot_map.end())
    {
        throw std::out_of_range(
            "ComponentHelper::GetSpotList ‑ residue is not supported");
    }
    return m_spot_map.at(residue);
}

const std::vector<Link> & ComponentHelper::GetLinkList(Residue residue)
{
    if (m_link_map.find(residue) == m_link_map.end())
    {
        throw std::out_of_range(
            "ComponentHelper::GetLinkList ‑ residue is not supported");
    }
    return m_link_map.at(residue);
}
