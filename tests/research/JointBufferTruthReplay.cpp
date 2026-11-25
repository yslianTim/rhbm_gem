#include "core/estimator/joint_component/LocalConditioning.hpp"
#include "core/estimator/joint_component/Problem.hpp"
#include "core/estimator/joint_component/StructuralPartition.hpp"
#include "data/io/detail/JointResultJson.hpp"

#include <boost/json.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <numbers>

namespace {
using namespace rhbm_gem;
namespace n = core::joint_component;
namespace j = boost::json;
namespace mp = boost::multiprecision;

j::value Finite(double x)
{
    return std::isfinite(x) ? j::value(x) : j::value();
}

// Mirror the replay's double kernel for Reference A. The native checker
// remains the source of measured gradient agreement and qualification.
n::BasisValues ReplayBasis(double square, double tau)
{
    constexpr double pi=std::numbers::pi;
    const double radius=std::sqrt(square), tau_squared=tau*tau;
    const double exponential=std::exp(-square/(2*tau_squared));
    const double gaussian=exponential/std::pow(2*pi*tau_squared,1.5);
    const double center=std::sqrt(2/(pi*tau_squared));
    return {gaussian,radius<1e-5 ? center : std::erf(radius/tau/std::sqrt(2.0))/radius,
        gaussian*(square/tau_squared-3),radius<1e-5 ? -center : -center*exponential};
}

// Research reconstruction exports rowwise quantities the native evidence
// reduces to scalars. Check those reductions against the actual replay.
std::pair<n::Vector,n::Vector> ReplayRows(const n::Domain & domain, n::VectorRef y,
    const n::Evaluation & primary, n::VectorRef tau, double scale)
{
    n::Vector prediction=n::Vector::Zero(y.size()), correction=prediction;
    const auto add_product=[&](Eigen::Index row, double a, double b) {
        const double product=a*b;
        for (const double term:{product,std::fma(a,b,-product)})
        {
            const double sum=prediction(row), next=sum+term;
            correction(row)+=std::abs(sum)>=std::abs(term) ? (sum-next)+term : (term-next)+sum;
            prediction(row)=next;
        }
    };
    for (Eigen::Index i=0;i<tau.size();++i)
        for (const auto & support:domain.atoms[static_cast<std::size_t>(i)])
        {
            const auto b=ReplayBasis(support.square,tau(i));
            add_product(support.row,primary.beta(2*i),b.gaussian);
            add_product(support.row,primary.beta(2*i+1),b.charge);
        }
    prediction+=correction;
    n::Vector gradients=n::Vector::Zero(tau.size());
    for (Eigen::Index i=0;i<tau.size();++i)
    {
        double compensation=0;
        for (const auto & support:domain.atoms[static_cast<std::size_t>(i)])
        {
            const auto b=ReplayBasis(support.square,tau(i));
            const double residual=prediction(support.row)-y(support.row);
            const double term=(primary.beta(2*i)*b.d_gaussian_d_log_tau+primary.beta(2*i+1)*b.d_charge_d_log_tau)*
                residual/scale/scale;
            const double adjusted=term-compensation, next=gradients(i)+adjusted;
            compensation=(next-gradients(i))-adjusted;
            gradients(i)=next;
        }
    }
    return {std::move(prediction),std::move(gradients)};
}

template <class Real> struct Basis
{
    Real g,k,dg,dk;
};

template <class Real> Basis<Real> ExactBasis(double squared, const Real & tau)
{
    const Real square(squared), radius=sqrt(square), variance=tau*tau;
    const Real pi=boost::math::constants::pi<Real>(), two_pi_variance=2*pi*variance;
    const Real exponential=exp(-square/(2*variance));
    const Real gaussian=exponential/(two_pi_variance*sqrt(two_pi_variance));
    const Real center=sqrt(Real(2)/pi)/tau;
    if (radius<Real("0.00001"))
        return {gaussian,center,gaussian*(square/variance-3),-center};
    const Real argument=radius/(sqrt(Real(2))*tau);
    return {gaussian,boost::math::erf(argument)/radius,gaussian*(square/variance-3),-center*exponential};
}

template <class Real> j::object Reference(const n::Domain & domain, n::VectorRef y,
    const n::Evaluation & primary, n::VectorRef tau, double scale, std::size_t target, int basis_source)
{
    std::vector<Real> residual(static_cast<std::size_t>(y.size()));
    std::vector<Basis<Real>> target_basis;
    for (std::size_t atom=0;atom<domain.atoms.size();++atom)
    {
        const auto i=static_cast<Eigen::Index>(atom);
        const Real width=basis_source==2 ? Real(exp(Real(primary.eta(i)))) : Real(tau(i));
        for (const auto & support:domain.atoms[atom])
        {
            Basis<Real> b;
            if (basis_source==2) b=ExactBasis(support.square,width);
            else
            {
                const auto values=basis_source==0 ? n::EvaluateKernel(support.square,tau(i),2.5) :
                    ReplayBasis(support.square,tau(i));
                b={Real(values.gaussian),Real(values.charge),Real(values.d_gaussian_d_log_tau),
                    Real(values.d_charge_d_log_tau)};
            }
            residual[static_cast<std::size_t>(support.row)]+=Real(primary.beta(2*i))*b.g+Real(primary.beta(2*i+1))*b.k;
            if (atom==target) target_basis.push_back(std::move(b));
        }
    }
    Real objective=0, prediction_difference=0;
    for (Eigen::Index row=0;row<y.size();++row)
    {
        auto & r=residual[static_cast<std::size_t>(row)];
        r-=Real(y(row));
        objective+=r*r;
        prediction_difference=std::max(prediction_difference,Real(abs(r-Real(primary.residual(row)))));
    }
    Real g=0, from_primary_residual=0, abs_terms=0, abs_gaussian=0, abs_charge=0, abs_derivative=0;
    std::size_t row_slot{};
    const auto target_column=2*static_cast<Eigen::Index>(target);
    const Real a(primary.beta(target_column)), c(primary.beta(target_column+1));
    for (const auto & support:domain.atoms[target])
    {
        const auto & b=target_basis[row_slot++];
        const Real u=a*b.dg, v=c*b.dk;
        const Real term=residual[static_cast<std::size_t>(support.row)]*(u+v);
        g+=term;
        from_primary_residual+=Real(primary.residual(support.row))*(u+v);
        abs_terms+=abs(term);
        abs_gaussian+=abs(u);
        abs_charge+=abs(v);
        abs_derivative+=abs(u+v);
    }
    const Real scale2=Real(scale)*scale;
    const auto decimal=[](const Real & x) {return x.str(0,std::ios_base::scientific);};
    return {{"gradient",decimal(g/scale2)}, {"gradient_using_primary_residual",decimal(from_primary_residual/scale2)},
        {"normalized_objective",decimal(objective/(2*scale2))},
        {"max_residual_difference_from_primary",decimal(prediction_difference)},
        {"gaussian_derivative_l1",decimal(abs_gaussian)}, {"charge_derivative_l1",decimal(abs_charge)},
        {"derivative_cancellation",abs_derivative==0 ? j::value() : j::value(decimal((abs_gaussian+abs_charge)/abs_derivative))},
        {"voxel_summation_cancellation",g==0 ? j::value() : j::value(decimal(abs_terms/abs(g)))}};
}
}

j::array DiagnoseTruthReplay(const core::JointProblem & problem, const core::JointFitResult & fit,
    const n::BcdExecutionPolicy & policy)
{
    const auto & data=core::JointProblemAccess::Get(problem);
    j::array result;
    for (std::size_t index=0;index<fit.components.size();++index)
    {
        const auto & c=fit.components[index];
        if (c.state || !c.candidate_diagnostic || !c.candidate_diagnostic->candidate_available) continue;
        const auto & candidate=*c.candidate_diagnostic;
        const auto & view=data.partition.components.at(index);
        JointParameterLayout layout;
        for (auto atom:data.layout.full_atoms)
            if (view.LocalAtom(static_cast<Eigen::Index>(atom))>=0) layout.full_atoms.push_back(atom);
        for (auto row:data.layout.informative_rows)
            if (view.LocalRow(static_cast<Eigen::Index>(row))>=0) layout.informative_rows.push_back(row);
        const auto prepared=n::MakeComponentProblemView(layout,view);
        const auto & domain=prepared.domain;
        auto context=data.context;
        context.atom_ids=prepared.atom_ids;
        context.row_ids=prepared.row_ids;
        context.width_bounds=policy.width_bounds;
        context.rank={prepared.parent_row_count,2*static_cast<Eigen::Index>(domain.atoms.size()),
            static_cast<Eigen::Index>(domain.atoms.size())};
        context.linear.rank_relative=context.rank.Relative(context.rank.design_columns);
        n::Vector y(domain.rows);
        for (Eigen::Index i=0;i<y.size();++i) y(i)=data.y(prepared.local_to_parent_rows[static_cast<std::size_t>(i)]);
        const n::VectorMap eta(candidate.candidate_eta.data(),static_cast<Eigen::Index>(candidate.candidate_eta.size()));
        const n::VectorMap beta(candidate.candidate_ac.data(),static_cast<Eigen::Index>(candidate.candidate_ac.size()));
        const n::Vector tau=eta.array().exp();
        const auto primary=n::EvaluateStateStreaming(domain,y,eta,beta,context);
        const auto replay=n::CheckIndependentReplay(domain,y,primary,context,nullptr,nullptr,n::PracticalJointConvergenceCriteria);
        const auto [independent_prediction,independent_gradient]=ReplayRows(domain,y,primary,tau,context.scale);
        const n::Vector independent_residual=independent_prediction-y;
        const double raw_gradient_difference=(primary.gradient-independent_gradient).lpNorm<Eigen::Infinity>();
        if (raw_gradient_difference!=replay.gradient_difference ||
            (primary.residual-independent_residual).lpNorm<Eigen::Infinity>()!=replay.residual_difference ||
            n::CalculateWidthBoundKkt(eta,independent_gradient,context.width_bounds).projected!=replay.independent_width_projected_kkt)
            throw std::logic_error("Rowwise diagnostics differ from native replay reductions.");
        j::object rec{{"trusted",false},{"component",c.id},{"stop",c.stop_reason},
            {"last_certification_sweep",c.endpoint_certification_sweep},{"last_search_sweep",c.search_sweeps},
            {"endpoint_reason",candidate.verification_reason},{"replay_reason",replay.reason},
            {"saved_normalized_objective",j::value_from(candidate.candidate_objective)},
            {"scale",context.scale},{"lower_eta",*policy.width_bounds.lower},{"upper_eta",*policy.width_bounds.upper},
            {"primary_ac_kkt",Finite(replay.primary_ac_kkt)},{"independent_ac_kkt",Finite(replay.independent_ac_kkt)},
            {"primary_projected_width_kkt",Finite(replay.primary_width_projected_kkt)},
            {"independent_projected_width_kkt",Finite(replay.independent_width_projected_kkt)},
            {"feasible",replay.feasible && replay.width_bound_feasible},
            {"prediction_difference",Finite(replay.prediction_difference)},
            {"residual_difference",Finite(replay.residual_difference)},
            {"objective_difference",Finite(replay.objective_difference)},
            {"primary_objective_unscaled",Finite(replay.primary_objective)},
            {"independent_objective_unscaled",Finite(replay.independent_objective)},
            {"prediction_cancellation",Finite(replay.cancellation_ratio)},
            {"comparison_flags",j::array{replay.prediction_passed,replay.residual_passed,replay.objective_passed,
                replay.kkt_agrees,replay.primary_width_kkt_passed,replay.independent_width_kkt_passed,replay.gradient_agrees}}};
        const double native_max_ratio_upper_bound=replay.gradient_difference/1e-13;
        rec["enclosure"]=j::object{{"atom",replay.maximum_gradient_difference_atom_id},
            {"primary",Finite(replay.maximum_gradient_difference_primary)},
            {"independent",Finite(replay.maximum_gradient_difference_independent)},
            {"difference",Finite(replay.gradient_difference)},{"bound",Finite(replay.maximum_gradient_difference_bound)},
            {"ratio_at_max_absolute_difference",Finite(replay.maximum_gradient_difference_ratio)},
            {"all_coordinate_ratio_upper_bound",Finite(native_max_ratio_upper_bound)},
            {"failed_atom",replay.gradient_enclosure_failure_atom_id}};
        const n::Vector projected=(eta-(eta-primary.gradient).cwiseMax(*policy.width_bounds.lower).
            cwiseMin(*policy.width_bounds.upper)).cwiseAbs();
        Eigen::Index limiting{};
        projected.maxCoeff(&limiting);
        j::array atoms, rows, blocks, references;
        for (Eigen::Index i=0;i<domain.rows;++i)
            rows.push_back(j::array{prepared.row_ids[static_cast<std::size_t>(i)],y(i),
                primary.residual(i)+y(i),primary.residual(i),independent_prediction(i),independent_residual(i)});
        for (Eigen::Index i=0;i<eta.size();++i)
        {
            const auto slot=static_cast<std::size_t>(i);
            const auto id=prepared.atom_ids[slot];
            j::array support;
            for (const auto & row:domain.atoms[slot]) support.push_back(j::array{row.row,row.square});
            atoms.push_back(j::object{{"id",id},{"eta",eta(i)},{"tau",tau(i)},{"A",beta(2*i)},{"C",beta(2*i+1)},
                {"amplitude_lower_face",beta(2*i)==0.0},{"primary_gradient",primary.gradient(i)},
                {"independent_gradient",independent_gradient(i)},{"primary_projected_kkt",projected(i)},
                {"width_lower_face",eta(i)-*policy.width_bounds.lower<=n::WidthBounds::RoundoffTolerance(*policy.width_bounds.lower)},
                {"width_upper_face",*policy.width_bounds.upper-eta(i)<=n::WidthBounds::RoundoffTolerance(*policy.width_bounds.upper)},
                {"support",std::move(support)}});
            if (i!=limiting && id!=replay.maximum_gradient_difference_atom_id) continue;
            j::object evidence{{"atom",id},{"limits_primary_width_kkt",i==limiting}};
            const auto compensation=n::AssessAcCompensation(domain,y,eta,beta,context,i);
            evidence["AC_conditioning"]=j::object{{"conditional_correlation",j::value_from(compensation.conditional_correlation)},
                {"cancellation_ratio",j::value_from(compensation.cancellation_ratio)},{"reason",compensation.reason}};
            evidence["reference_A_primary_basis"]=Reference<mp::cpp_dec_float_100>(domain,y,primary,tau,context.scale,slot,0);
            evidence["reference_A_independent_basis"]=Reference<mp::cpp_dec_float_100>(domain,y,primary,tau,context.scale,slot,1);
            evidence["reference_B_50"]=Reference<mp::cpp_dec_float_50>(domain,y,primary,tau,context.scale,slot,2);
            evidence["reference_B_100"]=Reference<mp::cpp_dec_float_100>(domain,y,primary,tau,context.scale,slot,2);
            references.push_back(std::move(evidence));
        }
        n::StructuralTopology topology(domain,context.atom_ids,prepared.local_layout);
        const auto partition=topology.BuildCores(n::ProductionBcdCoreAtoms);
        for (std::size_t i=0;i<partition.cores.size();++i)
        {
            j::array ids;
            for (auto atom:partition.cores[i].atoms) ids.emplace_back(prepared.atom_ids[static_cast<std::size_t>(atom)]);
            blocks.push_back(j::object{{"id",c.id+"/block-"+std::to_string(i+1)},{"core",std::move(ids)}});
        }
        rec["columns"]=j::array{"voxel","observation","primary_reconstructed_prediction","primary_residual",
            "independent_prediction","independent_residual"};
        rec["rows"]=std::move(rows);
        rec["atoms"]=std::move(atoms);
        rec["blocks"]=std::move(blocks);
        rec["coordinate_references"]=std::move(references);
        result.push_back(std::move(rec));
    }
    return result;
}
