// Included only by the offline benchmark; supplied-state rank, never profiling/search.
void RunRank(const n::Evaluation & state,const n::EvaluationContext & context,bool oracle,j::object & report,const char * output)
{
    n::Indices free;
    for(Eigen::Index k=0;k<state.beta.size();++k) if(k%2 || state.beta(k)>0) free.push_back(k);
    const auto p=static_cast<Eigen::Index>(free.size());
    if(oracle && p>1024) throw std::invalid_argument("Rank oracle is limited to 1024 free columns");
    n::Sparse design(state.x.rows(),p); std::vector<Eigen::Triplet<double>> entries;
    for(Eigen::Index col=0;col<p;++col)
    {
        const auto original=free[static_cast<std::size_t>(col)]; const double norm=state.x.col(original).norm();
        if(!(norm>0) || !std::isfinite(norm)) throw std::runtime_error("Invalid normalized free design");
        for(n::Sparse::InnerIterator v(state.x,original);v;++v) entries.emplace_back(v.row(),col,v.value()/norm);
    }
    design.setFromTriplets(entries.begin(),entries.end());
    report["free_columns"]=p; report["rank_rows"]=context.rank.rows; report["decision_columns"]=p;
    report["rank_backend"]=oracle ? "dense-oracle" : "spqr-bounds";
    report["not_run"]=j::array{"ac-solve","operator-search","reference","assessment","uncertainty"};
    report["stage"]="rank-factor"; Snapshot(output,report);
    const auto started=Clock::now();
    try {
        const auto factor=n::FreeDesignFactor::Fixed(design,free);
        report["rank_factor_seconds"]=Seconds(started);
        report["stage"]="rank-evaluation"; Snapshot(output,report);
        if(oracle)
        {
            const auto result=n::EvaluateRank(factor->Compact(),{context.rank,p});
            report["rank_result"]=SvdRecord(result);
        }
        else
        {
            const auto result=n::EvaluateFreeDesignRank(design,factor.get(),{context.rank,p});
            const auto number=[](double value)->j::value {return std::isfinite(value) ? j::value(value) : j::value(nullptr);};
            report["rank_result"]=j::object{
                {"status",result.status==n::FreeDesignRankStatus::FullRank ? "full-rank" : result.status==n::FreeDesignRankStatus::Deficient ? "deficient" : "unavailable"},
                {"reason",result.reason},{"rank_lower",result.rank_lower},{"rank_upper",result.rank_upper},
                {"exact_rank",result.status!=n::FreeDesignRankStatus::Unavailable && result.rank_lower==result.rank_upper ? j::value(result.rank_lower) : j::value(nullptr)},
                {"minimum_lower",number(result.minimum_lower)},{"maximum_lower",number(result.maximum_lower)},
                {"maximum_upper",number(result.maximum_upper)},{"threshold_lower",number(result.threshold_lower)},
                {"threshold_upper",number(result.threshold_upper)},{"reconstruction_error",number(result.reconstruction_error)},
                {"orthogonal_minimum",number(result.orthogonal_minimum)},{"witness_upper",number(result.witness_upper)},
                {"seconds",result.seconds},{"entries",result.entries},{"workspace_bytes",result.workspace_bytes}};
        }
    } catch(const std::runtime_error & error) {report["rank_result"]=j::object{{"status","unavailable"},{"reason",error.what()}};}
    report["rank_wall_seconds"]=Seconds(started);
    report["rank_compact_extractions"]=n::SparseWorkForTesting().compact_extractions;
    report["exported_factor_bytes"]=n::SparseWorkForTesting().factor_storage_bytes;
}
