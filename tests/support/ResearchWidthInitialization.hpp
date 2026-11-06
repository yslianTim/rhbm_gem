#pragma once
#include "core/detail/joint_component/WidthCalibration.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace width_research {
using namespace rhbm_gem;
enum class WidthSource { MDPDE, OLS };
inline const char * SourceName(WidthSource source)
{ return source==WidthSource::MDPDE ? "MDPDE" : "OLS"; }
inline bool ValidWidth(double tau) { return std::isfinite(tau) && tau>0.; }
inline const char * WidthValidity(double tau)
{ return !std::isfinite(tau) ? "nonfinite" : tau==0. ? "zero" : tau<0. ? "negative" : "valid"; }
inline bool ValidParameters(const std::array<double,3> & p)
{ return std::isfinite(p[0]) && p[0]>=0. && ValidWidth(p[1]) && std::isfinite(p[2]); }
inline bool CommonGate(const JointInitializationAtom & atom)
{
    return IsQualifiedFirstStageEstimate(atom.native_status,atom.first_stage_fit_qualified,
        atom.first_stage_qualification,atom.first_stage_recovery);
}
struct SelectedWidth
{
    std::string id;
    bool target{},width_valid{},parameters_valid{},common_gate{},calibration_donor{};
    double raw{},used{};
    std::string fallback_reason,donor_exclusion;
};
struct WidthSelection
{
    WidthSource source{};
    bool available{};
    std::vector<double> tau;
    std::vector<SelectedWidth> atoms;
    std::size_t fallback_donor_count{},fallback_count{};
    std::optional<double> fallback_median;
};

// Derived evidence only: never relabel OLS as MDPDE or change original records.
inline WidthSelection SelectWidths(const JointInitialization & original,WidthSource source)
{
    if(original.atoms.size()!=original.tau.size()) throw std::invalid_argument("Width source record size mismatch.");
    WidthSelection out; out.source=source; out.available=true;
    std::vector<double> fallback_donors;
    for(std::size_t i=0;i<original.atoms.size();++i)
    {
        const auto & atom=original.atoms[i];
        const auto & p=source==WidthSource::MDPDE ? atom.mdpde : atom.ols;
        const bool failure=atom.reason.starts_with("initialization-exception");
        if(atom.target && (source==WidthSource::MDPDE ? atom.reason=="valid-width" &&
            atom.original_tau && ValidWidth(*atom.original_tau) : !failure && ValidWidth(p[1])))
            fallback_donors.push_back(source==WidthSource::MDPDE ? *atom.original_tau : p[1]);
    }
    std::sort(fallback_donors.begin(),fallback_donors.end());
    out.fallback_donor_count=fallback_donors.size();
    if(!fallback_donors.empty())
    {
        const auto m=fallback_donors.size()/2;
        out.fallback_median=fallback_donors.size()%2 ? fallback_donors[m] :
            fallback_donors[m-1]+(fallback_donors[m]-fallback_donors[m-1])/2.;
    }
    for(std::size_t i=0;i<original.atoms.size();++i)
    {
        const auto & atom=original.atoms[i];
        const auto & p=source==WidthSource::MDPDE ? atom.mdpde : atom.ols;
        SelectedWidth selected;
        selected.id=atom.id; selected.target=atom.target; selected.raw=p[1];
        selected.width_valid=ValidWidth(p[1]); selected.parameters_valid=ValidParameters(p);
        selected.common_gate=CommonGate(atom);
        const bool failure=atom.reason.starts_with("initialization-exception");
        const bool production_fallback=source==WidthSource::MDPDE && atom.seed_source=="median-fallback";
        selected.used=source==WidthSource::MDPDE ? original.tau[i] : p[1];
        if(production_fallback) selected.fallback_reason="production-MDPDE-median-fallback";
        else if(source==WidthSource::OLS && atom.target && !selected.width_valid && out.fallback_median)
        { selected.used=*out.fallback_median; selected.fallback_reason="OLS-median-fallback"; }
        if(!ValidWidth(selected.used)) selected.fallback_reason="unavailable-no-valid-source-width";
        const bool fallback=production_fallback || selected.fallback_reason=="OLS-median-fallback";
        if(fallback) ++out.fallback_count;
        if(fallback) selected.donor_exclusion="median-fallback";
        else if(failure) selected.donor_exclusion="initialization-failure";
        else if(!selected.common_gate) selected.donor_exclusion="unqualified-first-stage";
        else if(!std::isfinite(p[1])) selected.donor_exclusion="non-finite-width";
        else if(!(p[1]>0.)) selected.donor_exclusion="non-positive-width";
        if(!atom.target) selected.donor_exclusion="halo";
        selected.calibration_donor=selected.donor_exclusion.empty();
        if(atom.target) out.available &= ValidWidth(selected.used);
        out.tau.push_back(selected.used); out.atoms.push_back(std::move(selected));
    }
    return out;
}

// A transient bridge used ONLY by the production projection helper. Its records
// retain original OLS/MDPDE evidence; selected Width metadata is labeled research.
inline JointInitialization ProjectionInput(const JointInitialization & original,const WidthSelection & selected)
{
    auto out=original; out.tau=selected.tau; out.valid=selected.available;
    out.data_scope="research-selected-widths-projection-only";
    for(std::size_t i=0;i<out.atoms.size();++i)
    {
        auto & atom=out.atoms[i]; const auto & width=selected.atoms[i];
        atom.original_tau=std::isfinite(width.raw) ? std::optional<double>(width.raw) : std::nullopt;
        atom.used_tau=ValidWidth(width.used) ? std::optional<double>(width.used) : std::nullopt;
        atom.seed_source=std::string("research-")+SourceName(selected.source)+":"+
            (width.fallback_reason.empty() ? "fitted" : width.fallback_reason);
    }
    return out;
}
}
