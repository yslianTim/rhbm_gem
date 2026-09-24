// Included by the offline benchmark inside its existing anonymous namespace.
// This bounded path never calls reference, assessment, or derivative reduction.
std::string fixed_mode,fixed_preconditioner;
std::filesystem::path fixed_state;
n::Vector FixedVector(const j::value & value)
{
    const auto & array=value.as_array(); n::Vector out(static_cast<Eigen::Index>(array.size()));
    for(std::size_t k=0;k<array.size();++k) out(static_cast<Eigen::Index>(k))=j::value_to<double>(array[k]);
    return out;
}
j::array FreeColumns(n::VectorRef beta)
{
    j::array out; for(Eigen::Index k=0;k<beta.size();++k) if(k%2 || beta(k)>0) out.push_back(k); return out;
}
j::object FixedWork()
{
    const auto & w=n::SparseWorkForTesting(); const auto & op=n::OperatorWorkForTesting();
    return {{"eigen_least_squares_seconds",w.least_squares_seconds},{"q_actions",w.q_actions},{"q_seconds",w.q_seconds},{"triangular_solves",w.triangular_solves},
        {"triangular_seconds",w.triangular_seconds},{"compact_extractions",w.compact_extractions},
        {"compact_seconds",w.compact_seconds},{"fixed_factorizations",w.fixed_factorizations},
        {"fixed_factor_seconds",w.fixed_factor_seconds},{"symbolic",w.symbolic},{"numeric",w.numeric},
        {"symbolic_seconds",w.symbolic_seconds},{"numeric_seconds",w.numeric_seconds},
        {"factor_r_nonzeros",w.factor_nonzeros},{"exported_factor_bytes",w.factor_storage_bytes ? j::value(w.factor_storage_bytes) : j::value(nullptr)},
        {"reference_solves",w.reference_solves},{"derivative_preparations",w.derivative_preparations},
        {"svd_seconds",w.free_design_svd_seconds},{"jacobi_retries",w.jacobi_retries},{"jacobi_retry_seconds",w.jacobi_retry_seconds},
        {"operator_preparation_seconds",op.preparation_seconds},{"operator_design_seconds",op.design_seconds},
        {"operator_factor_seconds",op.factor_seconds},{"operator_rank_seconds",op.rank_seconds},
        {"operator_compact_seconds",op.compact_seconds},{"operator_svd_seconds",op.svd_seconds},
        {"normal_actions",op.normals},{"normal_seconds",op.normal_seconds}};
}
void RunFixed(const n::Domain & domain,n::VectorRef y,const n::Vector & b,const n::EvaluationContext & context,const char * output)
{
    if(fixed_mode!="freeze" && fixed_mode!="composed" && fixed_mode!="normal") throw std::invalid_argument("Invalid fixed diagnostic mode");
    if(!fixed_preconditioner.empty() && fixed_preconditioner!="identity" && fixed_preconditioner!="diagonal" && fixed_preconditioner!="schwarz") throw std::invalid_argument("Invalid fixed preconditioner");
    if(b.size()>512) throw std::invalid_argument("Fixed diagnostics are bounded to 512 atoms");
    if(fixed_mode=="freeze")
    {
        const auto e=n::EvaluateProfile(domain,y,b.array().log(),false,&context);
        if(!e.valid || !e.certificate.kkt_passed || !n::CheckReplay(domain,y,e,context).passed) throw std::runtime_error("Cannot freeze uncertified state");
        Write(output,{{"stage","frozen"},{"eta",Values(e.eta)},{"beta",Values(e.beta)},
            {"free_columns",FreeColumns(e.beta)},{"scale",context.scale},{"rank_rows",context.rank.rows},
            {"design_columns",context.rank.design_columns},{"width_columns",context.rank.width_columns}});
        return;
    }
#ifdef PR4_BASELINE_DRIVER
    if(fixed_mode!="composed") throw std::invalid_argument("Baseline supports composed action only");
#endif
    const auto frozen=Read(fixed_state.c_str(),true).as_object();
    const auto eta=FixedVector(frozen.at("eta")),beta=FixedVector(frozen.at("beta"));
    if(eta.size()!=b.size() || beta.size()!=2*b.size() || frozen.at("free_columns")!=FreeColumns(beta) ||
        j::value_to<double>(frozen.at("scale"))!=context.scale || j::value_to<Eigen::Index>(frozen.at("rank_rows"))!=context.rank.rows ||
        j::value_to<Eigen::Index>(frozen.at("design_columns"))!=context.rank.design_columns ||
        j::value_to<Eigen::Index>(frozen.at("width_columns"))!=context.rank.width_columns) throw std::runtime_error("Frozen context mismatch");
    const auto basis_started=Clock::now(); const auto e=n::EvaluateState(domain,y,eta,beta,context);
    if(!e.valid) throw std::runtime_error("Invalid supplied state");
    j::object report{{"stage","operator"},{"mode",fixed_mode},{"basis_seconds",Seconds(basis_started)},
        {"state_control",j::object{{"eta",Values(e.eta)},{"beta",Values(e.beta)},{"free_columns",FreeColumns(e.beta)},
            {"scale",context.scale},{"rank_rows",context.rank.rows},{"residual",Values(e.residual/context.scale)},
            {"gradient",Values(e.gradient)},{"objective",.5*(e.residual/context.scale).squaredNorm()}}}};
    Snapshot(output,report);
    n::SparseWorkForTesting()={}; n::OperatorWorkForTesting()={}; n::SearchWorkForTesting()={};
    svd_records.clear(); n::CompactSvdCaptureForTesting()=Capture;
    {
    const n::ProfileJacobianOperator op(e,context);
    n::CompactSvdCaptureForTesting()={};
    report["rank"]=svd_records; report["valid"]=op.Valid(); report["reason"]=op.Reason();
    report["preparation_work"]=FixedWork();
    if(!op.Valid()) {report["stage"]="complete"; Snapshot(output,report); return;}
    const n::Vector v=n::Vector::LinSpaced(op.Columns(),-.3,.7);
    n::Vector w(op.Rows()); for(Eigen::Index r=0;r<w.size();++r) w(r)=std::sin(.17*static_cast<double>(r));
    report["apply"]=Values(op.Apply(v)); report["adjoint"]=Values(op.ApplyAdjoint(w));
    const auto normal=[&](n::VectorRef x)->n::Vector {
#ifndef PR4_BASELINE_DRIVER
        if(fixed_mode=="normal") return op.ApplyNormal(x);
#endif
        return op.ApplyAdjoint(op.Apply(x));
    };
    const auto q_before=n::SparseWorkForTesting().q_actions;
    const auto normal_start=Clock::now(); const n::Vector nv=normal(v);
    report["normal_action_seconds"]=Seconds(normal_start); report["normal_q_actions"]=n::SparseWorkForTesting().q_actions-q_before;
    report["normal"]=Values(nv);
    report["total_includes_gradient"]=true;
    report["timing_scope"]="independent preparation through prediction; audit actions and serialization excluded";
    report["operator_gradient"]=Values(op.ApplyAdjoint(e.residual/context.scale));
    } // Release the audit factor before independent timed steps.
    j::array steps;
    for(const auto * kind:{"identity","diagonal","schwarz"})
    {
        if(!fixed_preconditioner.empty() && fixed_preconditioner!=kind) continue;
        report["stage"]=kind; Snapshot(output,report);
        const auto wall_start=Clock::now();
        const auto preparation_start=Clock::now();
        const n::ProfileJacobianOperator step_op(e,context);
        const double preparation_seconds=Seconds(preparation_start);
        if(!step_op.Valid())
        {
            steps.push_back(j::object{{"kind",kind},{"valid",false},{"reason",step_op.Reason()},
                {"fixed_step_wall_seconds",Seconds(wall_start)}});
            report["steps"]=steps; continue;
        }
        const auto gradient_start=Clock::now();
        const n::Vector gradient=step_op.ApplyAdjoint(e.residual/context.scale);
        const double gradient_seconds=Seconds(gradient_start);
        const auto metric_start=Clock::now();
        const auto norms=n::WidthNorms(n::RawWidthDerivative(e),context.scale); const auto metric=n::WidthMetric(norms);
        const double metric_seconds=Seconds(metric_start);
        const n::PreconditionerContext pc{step_op.Identity(),n::PreconditionerSpace::Width,metric,1e-3};
        double partition_seconds{};
        std::unique_ptr<n::SchwarzModel> model; std::unique_ptr<n::SchwarzPreconditioner> inverse;
        const auto build_start=Clock::now();
        if(std::string(kind)=="schwarz")
        {
            const auto partition_start=Clock::now(); const auto partition=n::SearchPartition(domain,context);
            partition_seconds=Seconds(partition_start);
            model=std::make_unique<n::SchwarzModel>(*partition,e,context.scale,pc);
            inverse=std::make_unique<n::SchwarzPreconditioner>(*model,pc);
        }
        const double build_seconds=Seconds(build_start)-partition_seconds;
        const n::Vector diagonal=norms.array().square()+pc.damping*metric.array().square();
        double inverse_seconds{}; std::size_t inverse_actions{};
        const auto precondition=[&](n::VectorRef r)->n::Vector {
            n::WorkTimer timer(inverse_seconds); ++inverse_actions;
            if(inverse) return inverse->ApplyInverse(r,pc);
            if(std::string(kind)=="diagonal") return r.array()/diagonal.array();
            return r;
        };
        const auto action=[&](n::VectorRef x)->n::Vector {
            n::Vector result;
#ifndef PR4_BASELINE_DRIVER
            if(fixed_mode=="normal") result=step_op.ApplyNormal(x);
            else
#endif
                result=step_op.ApplyAdjoint(step_op.Apply(x));
            return result+(pc.damping*metric.array().square()*x.array()).matrix();
        };
        const auto solve_start=Clock::now(); auto step=n::SolvePcg(action,precondition,-gradient,metric);
        if(step.valid) step.predicted=-gradient.dot(step.step)-.5*step_op.Apply(step.step).squaredNorm();
        const double solve_seconds=Seconds(solve_start),wall_seconds=Seconds(wall_start);
        steps.push_back(j::object{{"kind",kind},{"valid",step.valid},{"reason",step.reason},{"step",Values(step.step)},
            {"iterations",step.iterations},{"true_residual",std::isfinite(step.relative_residual) ? j::value(step.relative_residual) : j::value(nullptr)},
            {"predicted",std::isfinite(step.predicted) ? j::value(step.predicted) : j::value(nullptr)},
            {"preparation_seconds",preparation_seconds},{"gradient_seconds",gradient_seconds},{"metric_seconds",metric_seconds},{"partition_seconds",partition_seconds},
            {"build_seconds",build_seconds},{"inverse_seconds",inverse_seconds},{"inverse_actions",inverse_actions},{"solve_seconds",solve_seconds},
            {"fixed_step_wall_seconds",wall_seconds},
            {"total_seconds",preparation_seconds+gradient_seconds+metric_seconds+partition_seconds+build_seconds+solve_seconds}});
        report["steps"]=steps;
    }
    report["factor_work"]=FixedWork(); report["search_work"]=SearchWork(); report["stage"]="complete"; Snapshot(output,report);
}
