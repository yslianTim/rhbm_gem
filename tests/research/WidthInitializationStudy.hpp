// Included by FirstStageSamplingStudy.cpp to reuse its serializers and endpoint audit.
j::object DescribeCalibration(const n::WidthBoundCalibration & c)
{
    j::object out{{"status",c.status},{"reason",c.unavailable_reason},{"identity",c.dataset_identity},
        {"donors",c.qualified_donor_count},{"coverage",c.qualified_donor_coverage},
        {"fallback",c.fallback_count},{"native",c.native_success_count},{"refined",c.refined_success_count},
        {"unqualified",c.unqualified_count},{"invalid_width",c.invalid_width_count},
        {"initial_violations",c.initial_violations4}};
    j::object exclusions; for(const auto & [reason,count]:c.excluded_donor_reasons) exclusions[reason]=count;
    out["exclusions"]=exclusions;
    if(c.tau_max4) out["bound"]=*c.tau_max4;
    if(c.p95) { out["p50"]=*c.median_width; out["p90"]=*c.p90; out["p95"]=*c.p95; out["p99"]=*c.p99; }
    return out;
}

j::object DescribeSelection(const width_research::WidthSelection & selected)
{
    j::array atoms;
    for(const auto & a:selected.atoms)
        atoms.push_back(j::object{{"id",a.id},{"target",a.target},{"raw",Number(a.raw)},{"used",Number(a.used)},
            {"width_valid",a.width_valid},{"width_validity",width_research::WidthValidity(a.raw)},
            {"parameters_valid",a.parameters_valid},{"common_MDPDE_qualification_gate",a.common_gate},
            {"calibration_donor",a.calibration_donor},{"donor_exclusion",a.donor_exclusion},
            {"fallback_reason",a.fallback_reason}});
    j::object out{{"source",width_research::SourceName(selected.source)},{"available",selected.available},
        {"tau",Values(selected.tau)},{"fallback_donors",selected.fallback_donor_count},
        {"fallback_count",selected.fallback_count},{"atoms",atoms}};
    if(selected.fallback_median) out["fallback_median"]=*selected.fallback_median;
    return out;
}

void RunWidthStudy(j::object & output,const JointInitialization & original,const JointProblem & problem,
    MapObject & map,ModelObject & model,const d::FittingWorkset & workset,const FitOptions & options,bool real)
{
    using width_research::WidthSource;
    const auto adapter_start=Clock::now();
    const std::array selections{width_research::SelectWidths(original,WidthSource::MDPDE),
        width_research::SelectWidths(original,WidthSource::OLS)};
    output["adapter_seconds"]=Seconds(adapter_start);
    const auto calibration_start=Clock::now();
    const std::array calibrations{width_research::CalibrateWidths(original,selections[0]),
        width_research::CalibrateWidths(original,selections[1])};
    output["calibration_seconds"]=Seconds(calibration_start);
    output["selections"]=j::object{{"MDPDE",DescribeSelection(selections[0])},{"OLS",DescribeSelection(selections[1])}};
    output["calibrations"]=j::object{{"MDPDE",DescribeCalibration(calibrations[0])},{"OLS",DescribeCalibration(calibrations[1])}};
    std::size_t common=0,m_only=0,o_only=0,neither=0;
    for(std::size_t i=0;i<original.atoms.size();++i)
    {
        if(!original.atoms[i].target) continue;
        const bool m=selections[0].atoms[i].calibration_donor,o=selections[1].atoms[i].calibration_donor;
        if(m && o) ++common; else if(m) ++m_only; else if(o) ++o_only; else ++neither;
    }
    output["donor_cohorts"]=j::object{{"common",common},{"MDPDE_only",m_only},{"OLS_only",o_only},{"neither",neither}};
    if(m_only || o_only)
        output["matched_cohort_calibrations"]=j::object{
            {"MDPDE",DescribeCalibration(width_research::CalibrateWidths(original,selections[0],true))},
            {"OLS",DescribeCalibration(width_research::CalibrateWidths(original,selections[1],true))}};
    double common_bound=0.;
    for(const auto & c:calibrations) if(c.status=="Available") common_bound=std::max(common_bound,*c.tau_max4);
    output["common_bound"]=common_bound;
    output["observation_scale"]=problem.ObservationScale();
    output["observations"]=Values(problem.Input().observations);
    j::array ids; for(const auto & id:problem.Input().atom_ids) ids.push_back(j::value(id));
    output["problem_atom_ids"]=ids;
    double initial_objective=0.,constant_objective=0.;
    for(double y:problem.Input().observations) initial_objective+=(y/problem.ObservationScale())*(y/problem.ObservationScale())/2.;
    for(auto row:JointProblemAccess::Get(problem).partition.constant_rows)
    {
        const double y=problem.Input().observations[static_cast<std::size_t>(row)]/problem.ObservationScale();
        constant_objective+=y*y/2.;
    }
    output["initial_objective"]=initial_objective;
    j::object runs,reference_fit;
    for(bool common_mode:{false,true})
    {
        // Synthetic common-bound control is needed only for crowded phase zero.
        if(common_mode && !real && output.at("scenario").as_string()!="crowded:0") continue;
        for(const std::string id:{"MM","OM","MO","OO"})
        {
            const std::size_t seed=id[0]=='M' ? 0 : 1,bound_source=id[1]=='M' ? 0 : 1;
            const auto & selected=selections[seed]; const auto & calibration=calibrations[bound_source];
            const auto name=id+(common_mode ? "-common" : "");
            j::object run{{"id",id},{"seed_source",width_research::SourceName(selected.source)},
                {"bound_source",width_research::SourceName(selections[bound_source].source)},
                {"mode",common_mode ? "common-research" : id=="MM" ? "production" : "source-research"},
                {"seed_available",selected.available}};
            if(!selected.available || calibration.status!="Available")
            {
                run["unavailable_reason"]=selected.available ? calibration.unavailable_reason : "seed-unavailable";
                runs[name]=run; continue;
            }
            auto bound=calibration;
            if(common_mode) { bound.tau_max4=common_bound; bound.eta_max4=std::log(common_bound); }
            const auto projection_start=Clock::now();
            const auto bounded=n::MakeBoundedInitialization(width_research::ProjectionInput(original,selected),bound);
            run["projection_seconds"]=Seconds(projection_start);
            run["bound"]=bounded.tau_max; run["initial_tau"]=Values(selected.tau);
            run["bounded_tau"]=Values(bounded.tau_values); run["projection_count"]=bounded.projection_count;
            run["maximum_projection"]=bounded.maximum_projection_magnitude;
            j::array projections;
            for(const auto & p:bounded.projections) if(p.projection_applied)
                projections.push_back(j::object{{"id",p.atom_id},{"raw",p.original_first_stage_tau ? Number(*p.original_first_stage_tau) : j::value(nullptr)},
                    {"used",p.original_initial_tau},{"bounded",p.bounded_initial_tau}});
            run["projections"]=projections;
            j::object first_sweeps; Clock::time_point solver_start=Clock::now();
            const n::JointProgressObserver observer=[&](const n::JointProgressEvent & event) {
                if(event.phase==n::JointProgressPhase::SolverConfigured) solver_start=Clock::now();
                if(event.phase==n::JointProgressPhase::SearchProgress && event.bcd && event.bcd->sweep==1)
                    first_sweeps[event.component_id]=j::object{{"objective",Number(event.bcd->objective)},
                        {"accepted_updates",event.bcd->local_accepted_updates},{"evaluations",event.bcd->profile_evaluations}};
            };
            const auto fit=[&] {
                if(id=="MM" && !common_mode)
                    return d::FitJointComponentsInitializedWithProgress(problem,original,observer);
                n::BcdExecutionPolicy policy; policy.eta_upper_bound=bound.eta_max4;
                return n::FitBcdComponents(problem,bounded.tau_values,policy,observer);
            }();
            run["solver_seconds"]=Seconds(solver_start);
            run["peak_rss_bytes"]=PeakRssBytes();
            run["joint"]=DescribeFit(fit); run["first_sweeps"]=first_sweeps;
            double first_objective=constant_objective;
            for(const auto & entry:first_sweeps) first_objective+=j::value_to<double>(entry.value().as_object().at("objective"));
            run["first_sweep_objective"]=first_objective;
            run["first_sweep_decrease"]=initial_objective-first_objective;
            int accepted=0,evaluations=0;
            for(const auto & c:fit.components) { accepted+=c.accepted_updates; evaluations+=c.profile_evaluations; }
            run["accepted_updates"]=accepted; run["search_evaluations"]=evaluations;
            if(fit.assembled_state)
            {
                if(reference_fit.empty()) reference_fit=run.at("joint").as_object();
                run["endpoint_audit"]=AuditEndpoint(problem,reference_fit,run.at("joint").as_object(),bounded.tau_max);
            }
            if(id=="MM" && !common_mode)
            {
                // Exact Production provenance permits the unchanged downstream contract.
                const auto downstream_start=Clock::now(); d::JointPostprocessingMetrics metrics;
                d::RunJointPostprocessing(map,model,problem,CaptureJointAnalysisResult(fit),options,metrics);
                std::size_t local=0,group=0;
                for(const auto * atom:workset.contributors)
                {
                    const auto view=AtomLocalPotentialView::For(*atom);
                    if(view.GetPostJointLocalResult() && view.GetPostJointLocalResult()->status==EvidenceStatus::Available) ++local;
                    if(view.GetGroupMemberResult()) ++group;
                }
                run["downstream"]=j::object{{"status","production-validated"},{"local_available",local},{"group_available",group},
                    {"peeling_samples",metrics.grid_peeling_sample_count},{"peeling_available",metrics.available_peeling_sample_count},
                    {"seconds",Seconds(downstream_start)}};
            }
            else run["downstream"]=j::object{{"status","not-published"},
                {"reason","research-source/common-bound provenance has no Production persistence contract"}};
            runs[name]=std::move(run);
        }
    }
    j::object paired_audits;
    for(const auto & pair:std::array<std::array<std::string,2>,4>{{
        {"OM","OO"},{"MM-common","OM-common"},{"MM","MM-common"},{"OM","OM-common"}}})
    {
        if(!runs.contains(pair[0]) || !runs.contains(pair[1])) continue;
        const auto & a=runs.at(pair[0]).as_object(), & b=runs.at(pair[1]).as_object();
        if(!a.contains("joint") || !b.contains("joint") || !a.at("joint").as_object().contains("tau") ||
            !b.at("joint").as_object().contains("tau")) continue;
        paired_audits[pair[0]+"-vs-"+pair[1]]=AuditEndpoint(problem,a.at("joint").as_object(),
            b.at("joint").as_object(),j::value_to<double>(b.at("bound")));
    }
    output["paired_endpoint_audits"]=std::move(paired_audits);
    output["runs"]=std::move(runs);
    output["qscore"]=CalculateQScoreForAtom(*model.GetAtomList().front(),map,model,1.,0.);
}
