#include "AminoAcidInfoHelper.hpp"
#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"

#include <iostream>
#include <stdexcept>

const std::unordered_map<Residue, std::vector<Element>> AminoAcidInfoHelper::m_element_map
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

    {Residue::CSX, {Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
                    Element::CARBON, Element::SULFUR, Element::OXYGEN}}
};

const std::unordered_map<Residue, std::vector<Remoteness>> AminoAcidInfoHelper::m_remoteness_map
{
    {Residue::ALA, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA}},

    {Residue::ARG, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::EPSILON,
                    Remoteness::ZETA, Remoteness::ETA, Remoteness::ETA}},

    {Residue::ASN, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::DELTA}},

    {Residue::ASP, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::DELTA}},

    {Residue::CYS, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA}},

    {Residue::GLN, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::EPSILON,
                    Remoteness::EPSILON}},

    {Residue::GLU, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::EPSILON,
                    Remoteness::EPSILON}},

    {Residue::GLY, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE}},

    {Residue::HIS, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::DELTA,
                    Remoteness::EPSILON, Remoteness::EPSILON}},

    {Residue::ILE, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::GAMMA, Remoteness::DELTA}},

    {Residue::LEU, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::DELTA}},

    {Residue::LYS, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::EPSILON,
                    Remoteness::ZETA}},

    {Residue::MET, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::EPSILON}},

    {Residue::PHE, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::DELTA,
                    Remoteness::EPSILON, Remoteness::EPSILON, Remoteness::ZETA}},

    {Residue::PRO, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA}},

    {Residue::SER, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA}},

    {Residue::THR, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::GAMMA}},

    {Residue::TRP, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::DELTA,
                    Remoteness::EPSILON, Remoteness::EPSILON, Remoteness::EPSILON, Remoteness::ZETA,
                    Remoteness::ZETA, Remoteness::ETA}},

    {Residue::TYR, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA, Remoteness::DELTA,
                    Remoteness::EPSILON, Remoteness::EPSILON, Remoteness::ZETA, Remoteness::ETA}},

    {Residue::VAL, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::GAMMA}},

    {Residue::CSX, {Remoteness::NONE, Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE,
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::DELTA}}
};

const std::unordered_map<Residue, std::vector<Branch>> AminoAcidInfoHelper::m_branch_map
{
    {Residue::ALA, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE}},

    {Residue::ARG, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::ONE, Branch::TWO}},

    {Residue::ASN, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::ONE, Branch::TWO}},

    {Residue::ASP, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::ONE, Branch::TWO}},

    {Residue::CYS, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE}},

    {Residue::GLN, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::NONE, Branch::ONE,
                    Branch::TWO}},

    {Residue::GLU, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::NONE, Branch::ONE,
                    Branch::TWO}},

    {Residue::GLY, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE}},

    {Residue::HIS, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::ONE, Branch::TWO,
                    Branch::ONE, Branch::TWO}},

    {Residue::ILE, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::ONE, Branch::TWO, Branch::ONE}},

    {Residue::LEU, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::ONE, Branch::TWO}},

    {Residue::LYS, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE}},

    {Residue::MET, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE}},

    {Residue::PHE, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::ONE, Branch::TWO,
                    Branch::ONE, Branch::TWO, Branch::NONE}},

    {Residue::PRO, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::NONE}},

    {Residue::SER, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE}},

    {Residue::THR, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::ONE, Branch::TWO}},

    {Residue::TRP, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::ONE, Branch::TWO,
                    Branch::ONE, Branch::TWO, Branch::THREE, Branch::TWO,
                    Branch::THREE, Branch::TWO}},

    {Residue::TYR, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::ONE, Branch::TWO,
                    Branch::ONE, Branch::TWO, Branch::NONE, Branch::NONE}},

    {Residue::VAL, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::ONE, Branch::TWO}},
    
    {Residue::CSX, {Branch::NONE, Branch::NONE, Branch::NONE, Branch::NONE,
                    Branch::NONE, Branch::NONE, Branch::NONE}}
};

const std::unordered_map<Residue, std::vector<double>> AminoAcidInfoHelper::m_amber95_partial_charge_map
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
    {Residue::VAL, { 0.597,-0.088,-0.416,-0.568, 0.299,-0.319,-0.319}},
    {Residue::CSX, { 0.597, 0.021,-0.416,-0.568,-0.123,-0.285, 0.000}}
};

const std::unordered_map<Residue, std::vector<double>> AminoAcidInfoHelper::m_buried_partial_charge_map
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
    {Residue::VAL, { 0.560, 0.215,-0.588,-0.661,-0.005,-0.362,-0.360}},
    {Residue::CSX, { 0.574, 0.186,-0.558,-0.694, 0.096,-0.285, 0.000}}
};

const std::unordered_map<Residue, std::vector<double>> AminoAcidInfoHelper::m_helix_partial_charge_map
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
    {Residue::VAL, { 0.559, 0.220,-0.602,-0.695,-0.013,-0.362,-0.362}},
    {Residue::CSX, { 0.571, 0.188,-0.577,-0.718, 0.086,-0.285, 0.000}}
};

const std::unordered_map<Residue, std::vector<double>> AminoAcidInfoHelper::m_sheet_partial_charge_map
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
    {Residue::VAL, { 0.563, 0.214,-0.587,-0.657, 0.000,-0.361,-0.356}},
    {Residue::CSX, { 0.584, 0.183,-0.550,-0.695, 0.103,-0.285, 0.000}}
};

size_t AminoAcidInfoHelper::GetAtomCount(Residue residue)
{
    return m_element_map.at(residue).size();
}

size_t AminoAcidInfoHelper::GetAtomCount(int residue)
{
    return m_element_map.at(static_cast<Residue>(residue)).size();
}

double AminoAcidInfoHelper::GetPartialCharge(
    Residue residue, Element element, Remoteness remoteness, Branch branch,
    Structure structure, bool verbose)
{
    using Key = std::uint32_t; // packed <Element, Remoteness, Branch>
    auto pack_key = [](Element e, Remoteness r, Branch b) noexcept -> Key
    {
        return  (static_cast<Key>(e) << 16) |
                (static_cast<Key>(r) <<  8) |
                 static_cast<Key>(b);
    };

    // one cache bucket per residue, lazily initialised on first use
    static std::unordered_map<Residue, std::unordered_map<Key, double>> cache;

    auto & residue_cache{ cache[residue] };
    try
    {
        if (residue_cache.empty()) // first request for this residue
        {
            auto atom_size{ m_element_map.at(residue).size() };
            const auto & element_list{ m_element_map.at(residue) };
            const auto & remoteness_list{ m_remoteness_map.at(residue) };
            const auto & branch_list{ m_branch_map.at(residue) };
            const auto & charge_list{ GetPartialChargeList(residue, structure) };

            // the four vectors should guaranteed aligned
            if (atom_size != element_list.size() ||
                atom_size != remoteness_list.size() ||
                atom_size != branch_list.size() ||
                atom_size != charge_list.size())
            {
                throw std::range_error(
                    "AminoAcidInfoHelper::GetPartialCharge ‑ the four vectors are not aligned");
            }

            residue_cache.reserve(atom_size);
            for (std::size_t i = 0; i < atom_size; i++)
            {
                residue_cache.emplace(
                    pack_key(element_list[i], remoteness_list[i], branch_list[i]), charge_list[i]);
            }
        }
    }
    catch(const std::exception & except)
    {
        if (verbose == true) std::cout << except.what() << std::endl;
        return 0.0;
    }
    
    const Key key{ pack_key(element, remoteness, branch) };
    if (residue_cache.find(key) != residue_cache.end())
    {
        return residue_cache.at(key);
    }

    if (verbose == true) std::cout << "No partial charge data for this atom." << std::endl;
    return 0.0;
}

const std::vector<double> & AminoAcidInfoHelper::GetPartialChargeList(
    Residue residue, Structure structure)
{
    if (structure == Structure::FREE)
    {
        return m_buried_partial_charge_map.at(residue);
        //return m_amber95_partial_charge_map.at(residue); // Amber95 Table Test
    }
    else if (structure == Structure::SHEET)
    {
        return m_sheet_partial_charge_map.at(residue);
        //return m_amber95_partial_charge_map.at(residue); // Amber95 Table Test
    }
    else if (structure >= Structure::HELX_P && structure < Structure::TURN_P)
    {
        return m_helix_partial_charge_map.at(residue);
        //return m_amber95_partial_charge_map.at(residue); // Amber95 Table Test
    }
    else
    {
        throw std::out_of_range(
            "AminoAcidInfoHelper::GetPartialChargeList ‑ structure is not supported");
    }
}