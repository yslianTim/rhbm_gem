#include "AminoAcidInfoHelper.hpp"
#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"

#include <iostream>
#include <stdexcept>

const std::unordered_map<Residue, size_t> AminoAcidInfoHelper::m_atom_count_map
{
    {Residue::ALA,  5}, {Residue::ARG, 11}, {Residue::ASN,  8}, {Residue::ASP,  8},
    {Residue::CYS,  6}, {Residue::GLN,  9}, {Residue::GLU,  9}, {Residue::GLY,  4},
    {Residue::HIS, 10}, {Residue::ILE,  8}, {Residue::LEU,  8}, {Residue::LYS,  9},
    {Residue::MET,  8}, {Residue::PHE, 11}, {Residue::PRO,  7}, {Residue::SER,  6},
    {Residue::THR,  7}, {Residue::TRP, 14}, {Residue::TYR, 12}, {Residue::VAL,  7}
};

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
                    Element::CARBON, Element::CARBON, Element::CARBON}}
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
                    Remoteness::BETA, Remoteness::GAMMA, Remoteness::GAMMA}}
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
                    Branch::NONE, Branch::ONE, Branch::TWO}}
};

const std::unordered_map<Residue, std::vector<double>> AminoAcidInfoHelper::m_buried_partial_charge_map
{
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

size_t AminoAcidInfoHelper::GetAtomCount(Residue residue)
{
    return m_atom_count_map.at(residue);
}

size_t AminoAcidInfoHelper::GetAtomCount(int residue)
{
    return m_atom_count_map.at(static_cast<Residue>(residue));
}

double AminoAcidInfoHelper::GetPartialCharge(
    Residue residue, Element element, Remoteness remoteness, Branch branch, bool verbose)
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
            auto atom_size{ m_atom_count_map.at(residue) };
            const auto & element_list{ m_element_map.at(residue) };
            const auto & remoteness_list{ m_remoteness_map.at(residue) };
            const auto & branch_list{ m_branch_map.at(residue) };
            const auto & charge_list{ m_buried_partial_charge_map.at(residue) };

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
    const auto it{ residue_cache.find(key) };
    if (it != residue_cache.end())
    {
        return it->second;
    }

    return 0.0;
}