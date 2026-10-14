#include "Problem.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace rhbm_gem::core::joint_component {
JointParameterLayout BuildParameterLayout(const JointProblemInput & input)
{
    JointParameterLayout out;
    out.full_atoms.resize(input.atom_ids.size());
    std::iota(out.full_atoms.begin(),out.full_atoms.end(),0);
    out.informative_rows.resize(input.observations.size());
    std::iota(out.informative_rows.begin(),out.informative_rows.end(),0);
    return out;
}
std::vector<HaloSupportDiagnostic> DiagnoseHaloSupport(const JointProblemInput & input,
    const JointParameterLayout & layout,const std::vector<double> & widths,
    const std::vector<Eigen::Index> * component_atoms,Eigen::Index rank_context_rows)
{
    std::vector<HaloSupportDiagnostic> out;
    if(!input.selection_domain) return out;
    const auto is_in_component=[&](std::size_t atom) {
        if(!component_atoms) return true;
        return std::binary_search(component_atoms->begin(),component_atoms->end(),static_cast<Eigen::Index>(atom));
    };
    for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        if(std::binary_search(input.selection_domain->target_indices.begin(),
            input.selection_domain->target_indices.end(),atom) || !is_in_component(atom)) continue;
        HaloSupportDiagnostic diagnostic; diagnostic.atom_index=atom; diagnostic.atom_id=input.atom_ids[atom];
        const auto & support=input.support[atom];
        diagnostic.original_observation_support_count=support.size();
        diagnostic.singleton=support.size()==1;
        diagnostic.insufficient_full_abc_dimension=support.size()<3;
        std::vector<double> distances; distances.reserve(support.size());
        for(const auto & entry:support)
        {
            diagnostic.profiled_observation_support_count+=std::binary_search(
                layout.informative_rows.begin(),layout.informative_rows.end(),entry.row);
            distances.push_back(entry.squared_distance);
        }
        std::sort(distances.begin(),distances.end());
        distances.erase(std::unique(distances.begin(),distances.end()),distances.end());
        diagnostic.distinct_squared_distance_count=distances.size();
        diagnostic.exact_common_response=!support.empty() && distances.size()==1;
        const auto nonzero_columns=std::min<std::size_t>(support.size(),2);
        diagnostic.ac_structural_rank_upper_bound=nonzero_columns;
        if(atom<widths.size() && widths[atom]>0 && std::isfinite(widths[atom]) && !support.empty())
        {
            Eigen::MatrixXd design(static_cast<Eigen::Index>(support.size()),2);
            for(std::size_t row=0;row<support.size();++row)
            {
                const auto basis=EvaluateKernel(support[row].squared_distance,widths[atom],2.5);
                design(static_cast<Eigen::Index>(row),0)=basis.gaussian;
                design(static_cast<Eigen::Index>(row),1)=basis.charge;
            }
            const auto spectrum=ComputeSpectrum(design,
                {rank_context_rows>=0 ? rank_context_rows : static_cast<Eigen::Index>(input.observations.size()),2,0},2,true);
            if(spectrum.available)
            {
                diagnostic.gaussian_column_nonzero_norm=spectrum.column_norms(0)>0;
                diagnostic.charge_column_nonzero_norm=spectrum.column_norms(1)>0;
                diagnostic.fixed_width_ac_numerical_rank=static_cast<std::size_t>(spectrum.rank);
                diagnostic.column_normalized_condition_indicator=spectrum.condition;
                diagnostic.ac_structural_rank_upper_bound=std::min<std::size_t>(support.size(),
                    static_cast<std::size_t>(*diagnostic.gaussian_column_nonzero_norm)+
                    static_cast<std::size_t>(*diagnostic.charge_column_nonzero_norm));
            }
        }
        out.push_back(std::move(diagnostic));
    }
    return out;
}
}
