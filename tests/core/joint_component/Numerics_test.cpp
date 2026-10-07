#include <gtest/gtest.h>
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/joint_component/SparseFactor.hpp"
#include "core/detail/joint_component/CompactSvd.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointOperatorWorkload.hpp"
#include "support/JointTestNumerics.hpp"
#include <unsupported/Eigen/NonLinearOptimization>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace {
namespace m=second_stage_test::matched;
namespace p=m::joint_abc;
namespace n=p::runtime;
using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;
struct SvdMode
{
    p::runtime::CompactSvdMode saved{p::runtime::CompactSvdModeForTesting()};
    explicit SvdMode(p::runtime::CompactSvdMode mode) {p::runtime::CompactSvdModeForTesting()=mode;}
    ~SvdMode() {p::runtime::CompactSvdModeForTesting()=saved;}
};
struct Sample
{
    m::unique_grid::Grid grid;
    std::vector<m::Atom> atoms{{{-.6,0,0},2,.42,-.3},{{.7,.2,0},1.5,.67,.4}};
    Vector y;
    Sample()
    {
        for (int x=-6;x<=6;++x) for (int z=-3;z<=3;++z) for (int y0=-4;y0<=4;++y0)
        {
            m::unique_grid::Voxel voxel; voxel.position={x*.4,y0*.4,z*.4}; voxel.index=grid.voxels.size();
            grid.voxels.push_back(voxel);
        }
        y.resize(static_cast<Eigen::Index>(grid.voxels.size()));
        for (Eigen::Index k=0;k<y.size();++k) y(k)=m::unique_grid::Direct(grid.voxels[static_cast<std::size_t>(k)].position,atoms,2.5);
    }
};
Vector Read(const boost::json::value & v)
{
    Vector out(static_cast<Eigen::Index>(v.as_array().size()));
    for (Eigen::Index k=0;k<out.size();++k) out(k)=boost::json::value_to<double>(v.at(static_cast<std::size_t>(k)));
    return out;
}
// Independent small direct joint solve, with numerical derivatives of A/C/log-B.
// This comparison uses an interior solution; boundary behavior is tested separately.
struct Direct
{
    const Sample & sample;
    int values() const {return static_cast<int>(sample.y.size());}
    int operator()(const Vector & x,Vector & residual) const
    {
        auto atoms=sample.atoms;
        for (std::size_t a=0;a<atoms.size();++a)
        {atoms[a].amplitude=x(static_cast<Eigen::Index>(3*a)); atoms[a].charge=x(static_cast<Eigen::Index>(3*a+1)); atoms[a].width=std::exp(x(static_cast<Eigen::Index>(3*a+2)));}
        residual.resize(sample.y.size());
        for (Eigen::Index k=0;k<residual.size();++k)
            residual(k)=m::unique_grid::Direct(sample.grid.voxels[static_cast<std::size_t>(k)].position,atoms,2.5)-sample.y(k);
        return 0;
    }
    int df(const Vector & x,Matrix & jacobian) const
    {
        jacobian.resize(sample.y.size(),x.size());
        for (Eigen::Index k=0;k<x.size();++k)
        {
            Vector plus=x,minus=x,rp,rm; plus(k)+=1e-5; minus(k)-=1e-5;
            (*this)(plus,rp); (*this)(minus,rm); jacobian.col(k)=(rp-rm)/2e-5;
        }
        return 0;
    }
};

struct CompactJacobianDiagnostics
{
    double ztp_relative{},ztr_relative{},ztr_relative_all{},gram_absolute{},gram_relative{},gram_elementwise{};
    double orthogonal_gram_relative{},spectrum_relative{},weak_projector_difference{};
    double response_norm_difference{},response_gradient_difference{},correction_difference{};
    double ztp_absolute{},ztr_absolute{},ztr_scaled{};
    Eigen::Index exact_rank{},compact_rank{};
    bool projection_normalized{},stationarity_normalized{},rank_parity{},candidate_compared{};
};
struct CompactJacobianSummary
{
    CompactJacobianDiagnostics maximum;
    bool rank_parity{true};
};
void CheckCompactJacobianIdentity(const n::Evaluation & e,const n::EvaluationContext & context,
    const std::string & name,bool require_orthogonality,bool compare_compact_candidate,bool require_candidate_parity,
    CompactJacobianSummary & summary)
{
    SCOPED_TRACE(name);
    const auto d=n::PrepareDerivative(e,context.scale,&context);
    ASSERT_TRUE(d.valid)<<d.reason;
    const auto reduced=n::ReduceDerivative(d,e.residual,true);
    ASSERT_TRUE(reduced.valid)<<reduced.reason;

    const auto m=e.eta.size(),p=d.free_design.cols(),rows=e.residual.size();
    n::Matrix cross=n::Matrix::Zero(p,m); n::Vector ztr=n::Vector::Zero(p),projected_rhs=n::Vector::Zero(m);
    std::optional<n::TiledQR> design_qr;
    if(compare_compact_candidate) design_qr.emplace(p,1,"identity-design");
    n::Matrix projected,jacobian;
    for(Eigen::Index first=0;first<rows;first+=n::derivative_tile_rows)
    {
        const auto count=std::min(n::derivative_tile_rows,rows-first);
        d.Rows(first,count,projected,jacobian);
        const auto design=d.free_design.middleRows(first,count);
        cross.noalias()+=design.transpose()*projected;
        ztr.noalias()+=design.transpose()*e.residual.segment(first,count);
        if(compare_compact_candidate)
        {
            n::Matrix rhs(count,1); rhs.col(0)=e.residual.segment(first,count)/d.scale;
            projected_rhs.noalias()+=projected.transpose()*rhs.col(0);
            design_qr->Append(n::Matrix(design),rhs);
        }
    }

    CompactJacobianDiagnostics result;
    result.ztp_absolute=cross.norm();
    const double z_norm=d.free_design.norm(),p_norm=reduced.projected.norm(),residual_norm=e.residual.norm();
    const double orthogonal_scale=z_norm*p_norm;
    result.projection_normalized=orthogonal_scale>1e-30;
    if(result.projection_normalized) result.ztp_relative=result.ztp_absolute/orthogonal_scale;
    const double stationarity_scale=z_norm*residual_norm;
    result.ztr_absolute=ztr.norm();
    result.stationarity_normalized=residual_norm>1e-12*context.scale && stationarity_scale>1e-30;
    if(stationarity_scale>1e-30) result.ztr_relative=result.ztr_absolute/stationarity_scale;
    result.ztr_relative_all=result.ztr_relative;
    result.ztr_scaled=z_norm*context.scale>1e-30 ? result.ztr_absolute/(z_norm*context.scale) : result.ztr_absolute;

    const n::Matrix ztz=n::Matrix(d.free_design.transpose()*d.free_design);
    const n::Matrix correction=d.correction/d.scale;
    const n::Matrix projected_gram=reduced.projected.transpose()*reduced.projected;
    const n::Matrix exact_gram=reduced.jacobian.transpose()*reduced.jacobian;
    const n::Matrix full_gram=projected_gram-cross.transpose()*correction-correction.transpose()*cross+
        correction.transpose()*ztz*correction;
    const n::Matrix orthogonal_gram=projected_gram+correction.transpose()*ztz*correction;
    const n::Matrix gram_error=full_gram-exact_gram;
    result.gram_absolute=gram_error.norm();
    result.gram_relative=result.gram_absolute/std::max(1e-30,exact_gram.norm());
    result.gram_elementwise=gram_error.cwiseAbs().maxCoeff();
    result.orthogonal_gram_relative=(orthogonal_gram-exact_gram).norm()/std::max(1e-30,exact_gram.norm());
    EXPECT_LE(result.gram_relative,1e-8);
    if(require_orthogonality)
    {
        ASSERT_TRUE(result.projection_normalized);
        EXPECT_LE(result.ztp_relative,1e-8);
        EXPECT_LE(result.ztr_scaled,1e-8);
    }

    const n::RankRequest request{{context.rank.rows,0,0},m};
    const auto exact=n::EvaluateRank(reduced.jacobian,request,nullptr,n::CompactSvdVectors::Right);
    ASSERT_TRUE(exact.valid);
    result.exact_rank=exact.rank;
    if(compare_compact_candidate)
    {
        const n::Matrix & projected_factor=reduced.projected;
        const n::Vector projected_target=projected_factor.transpose().triangularView<Eigen::Lower>().solve(projected_rhs);
        n::Matrix compact_design(projected_factor.rows()+design_qr->r.rows(),m);
        compact_design.topRows(projected_factor.rows())=projected_factor;
        compact_design.bottomRows(design_qr->r.rows())=-(design_qr->r*correction);
        n::Matrix compact_response(compact_design.rows(),1);
        compact_response.topRows(projected_target.size())=projected_target;
        compact_response.bottomRows(design_qr->target.rows())=design_qr->target;
        n::TiledQR compact_qr(m,1,"identity-orthogonal-stack");
        compact_qr.Append(compact_design,compact_response);
        const auto compact=n::EvaluateRank(compact_qr.r,request,nullptr,n::CompactSvdVectors::Right);
        ASSERT_TRUE(compact.valid);
        result.candidate_compared=true;
        result.compact_rank=compact.rank;
        result.rank_parity=exact.rank==compact.rank;
        result.spectrum_relative=(exact.singular_values-compact.singular_values).lpNorm<Eigen::Infinity>()/
            std::max(1e-30,exact.singular_values(0));
        const Eigen::Index weak=std::min<Eigen::Index>(3,exact.right_vectors.cols());
        const n::Matrix exact_weak=exact.right_vectors.rightCols(weak),compact_weak=compact.right_vectors.rightCols(weak);
        result.weak_projector_difference=(exact_weak*exact_weak.transpose()-compact_weak*compact_weak.transpose()).norm()/
            std::sqrt(static_cast<double>(weak));

        const n::Vector exact_rhs=-reduced.response,compact_rhs=-compact_qr.target.col(0);
        const auto exact_solve=n::EvaluateRank(reduced.jacobian,request,&exact_rhs);
        const auto compact_solve=n::EvaluateRank(compact_qr.r,request,&compact_rhs);
        ASSERT_TRUE(exact_solve.valid && compact_solve.valid);
        result.correction_difference=(exact_solve.solution-compact_solve.solution).norm()/
            std::max(1.0,exact_solve.solution.norm());
        const n::Vector exact_gradient=reduced.jacobian.transpose()*reduced.response;
        const n::Vector compact_gradient=compact_qr.r.transpose()*compact_qr.target.col(0);
        result.response_gradient_difference=(exact_gradient-compact_gradient).norm()/std::max(1.0,exact_gradient.norm());
        result.response_norm_difference=std::abs(reduced.response.norm()-compact_qr.target.norm())/
            std::max(1.0,reduced.response.norm());
        if(require_candidate_parity)
        {
            EXPECT_TRUE(result.rank_parity);
            EXPECT_LE(result.spectrum_relative,1e-8);
            EXPECT_LE(result.response_gradient_difference,1e-8);
            EXPECT_LE(result.correction_difference,1e-7);
        }
        summary.rank_parity &= result.rank_parity;
    }
    auto & maximum=summary.maximum;
    maximum.ztp_relative=std::max(maximum.ztp_relative,result.ztp_relative);
    if(result.stationarity_normalized) maximum.ztr_relative=std::max(maximum.ztr_relative,result.ztr_relative);
    maximum.ztr_relative_all=std::max(maximum.ztr_relative_all,result.ztr_relative_all);
    maximum.ztr_scaled=std::max(maximum.ztr_scaled,result.ztr_scaled);
    maximum.ztp_absolute=std::max(maximum.ztp_absolute,result.ztp_absolute);
    maximum.gram_absolute=std::max(maximum.gram_absolute,result.gram_absolute);
    maximum.gram_relative=std::max(maximum.gram_relative,result.gram_relative);
    maximum.gram_elementwise=std::max(maximum.gram_elementwise,result.gram_elementwise);
    if(result.candidate_compared)
    {
        maximum.orthogonal_gram_relative=std::max(maximum.orthogonal_gram_relative,result.orthogonal_gram_relative);
        maximum.spectrum_relative=std::max(maximum.spectrum_relative,result.spectrum_relative);
        maximum.weak_projector_difference=std::max(maximum.weak_projector_difference,result.weak_projector_difference);
        maximum.response_norm_difference=std::max(maximum.response_norm_difference,result.response_norm_difference);
        maximum.response_gradient_difference=std::max(maximum.response_gradient_difference,result.response_gradient_difference);
        maximum.correction_difference=std::max(maximum.correction_difference,result.correction_difference);
    }
}
}

TEST(JointTestNumericsTest, LogWidthBasisIncludesCenterNearZeroAndFixedCutoff)
{
    for (double r:{0.,.5e-5,1e-5,1.1e-5,.4,2.5,2.50000001})
    {
        const double b=.53,h=1e-5;
        const auto value=m::EvaluateBasis(r*r,b,2.5),plus=m::EvaluateBasis(r*r,b*std::exp(h),2.5),minus=m::EvaluateBasis(r*r,b*std::exp(-h),2.5);
        EXPECT_NEAR(value.gaussian_log_width,(plus.gaussian-minus.gaussian)/(2*h),1e-8);
        EXPECT_NEAR(value.charge_log_width,(plus.charge-minus.charge)/(2*h),1e-8);
        if (r>2.5) EXPECT_EQ(value.gaussian+value.charge+value.gaussian_log_width+value.charge_log_width,0);
    }
}

TEST(JointTestNumericsTest, FullJacobianAndEnvelopeGradientAtNonzeroResidual)
{
    Sample sample; const p::Domain domain(sample.grid,sample.atoms);
    for (Eigen::Index k=0;k<sample.y.size();++k) sample.y(k)+=.03*std::sin(static_cast<double>(k));
    const Vector eta=Eigen::Vector2d(.55,.51).array().log();
    const auto e=p::Evaluate(domain,sample.y,eta); ASSERT_TRUE(e.valid);
    const double scale=std::max(1.0,sample.y.norm()); const auto d=p::Differentiate(e,scale);
    ASSERT_TRUE(d.valid); EXPECT_GT((d.jacobian-d.projected).norm(),1e-3);
    for (Eigen::Index k=0;k<2;++k)
    {
        Vector plus=eta,minus=eta; plus(k)+=1e-5; minus(k)-=1e-5;
        const auto ep=p::Evaluate(domain,sample.y,plus,true),em=p::Evaluate(domain,sample.y,minus,true);
        ASSERT_TRUE(ep.valid && em.valid);
        const Vector finite=(ep.residual-em.residual)/(2e-5*scale);
        EXPECT_LT((finite-d.jacobian.col(k)).norm()/finite.norm(),1e-7);
        EXPECT_NEAR(e.gradient(k),(ep.residual.squaredNorm()-em.residual.squaredNorm())/(4e-5*scale*scale),1e-9);
    }
}

TEST(JointTestNumericsTest, RecoveryAgreesWithDirectJointAndFixedB)
{
    Sample sample; const p::Domain domain(sample.grid,sample.atoms);
    const Vector widths=Eigen::Vector2d(.48,.59);
    const auto fit=p::Fit(domain,sample.y,widths);
    ASSERT_EQ(fit.at("runtime_convergence"),"passed")<<boost::json::serialize(fit);
    const auto beta=Read(fit.at("primary").at("beta")),b=Read(fit.at("primary").at("b"));
    EXPECT_LT((beta-Eigen::Vector4d(2,-.3,1.5,.4)).norm(),1e-10);
    EXPECT_LT((b-Eigen::Vector2d(.42,.67)).norm(),1e-10);
    EXPECT_EQ(fit.at("directional_evaluations"),0);
    const auto e=p::Evaluate(domain,sample.y,widths.array().log());
    const auto fixed=p::Evaluate(domain,sample.y,widths.array().log(),true);
    ASSERT_TRUE(fixed.valid); EXPECT_LT((e.beta-fixed.beta).norm(),1e-10);
    Direct direct{sample}; Eigen::LevenbergMarquardt<Direct> lm(direct);
    Vector joint(6); joint<<e.beta(0),e.beta(1),std::log(widths(0)),e.beta(2),e.beta(3),std::log(widths(1));
    lm.parameters.ftol=1e-14; lm.parameters.xtol=1e-12; lm.parameters.gtol=1e-12; lm.minimize(joint);
    EXPECT_NEAR(joint(0),beta(0),1e-10); EXPECT_NEAR(joint(3),beta(2),1e-10);
    EXPECT_NEAR(std::exp(joint(2)),b(0),1e-10); EXPECT_NEAR(std::exp(joint(5)),b(1),1e-10);
    EXPECT_GT(std::abs(b(0)-b(1)),.1); // Separate widths, not a shared B.
}

TEST(JointTestNumericsTest, BoundaryAndChangingActiveFaceRebuildDerivative)
{
    Sample sample; sample.atoms.resize(1); sample.atoms[0].amplitude=.05; sample.atoms[0].charge=-1;
    for (Eigen::Index k=0;k<sample.y.size();++k) sample.y(k)=m::unique_grid::Direct(sample.grid.voxels[static_cast<std::size_t>(k)].position,sample.atoms,2.5);
    const p::Domain domain(sample.grid,sample.atoms); bool active=false,interior=false;
    for (double width:{.25,.42,.8})
    {
        const Vector eta=Vector::Constant(1,std::log(width)); const auto e=p::Evaluate(domain,sample.y,eta);
        ASSERT_TRUE(e.valid); active |= e.beta(0)==0; interior |= e.beta(0)>0;
        const double scale=std::max(1.0,sample.y.norm()); const auto d=p::Differentiate(e,scale); ASSERT_TRUE(d.valid);
        const auto plus=p::Evaluate(domain,sample.y,eta.array()+1e-5,true),minus=p::Evaluate(domain,sample.y,eta.array()-1e-5,true);
        ASSERT_EQ(plus.certificate.at("active_atoms"),minus.certificate.at("active_atoms"));
        EXPECT_LT(((plus.residual-minus.residual)/(2e-5*scale)-d.jacobian.col(0)).norm(),1e-8);
    }
    EXPECT_TRUE(active && interior);
}

TEST(JointTestNumericsTest, RankFailureInvalidBAndUninformativeWidthCannotQualify)
{
    Sample sample; auto duplicate=sample.atoms; duplicate[1].position=duplicate[0].position;
    const p::Domain bad(sample.grid,duplicate),good(sample.grid,sample.atoms);
    const Vector eta=Vector::Constant(2,std::log(.5));
    const auto invalid=p::Evaluate(bad,sample.y,eta); EXPECT_FALSE(invalid.valid);
    EXPECT_FALSE(p::Differentiate(invalid,1).valid);
    EXPECT_NE(p::Fit(bad,sample.y,Vector::Constant(2,.5)).at("runtime_convergence"),"passed");
    EXPECT_FALSE(p::Evaluate(good,sample.y,Vector::Constant(2,1000)).valid);
    EXPECT_FALSE(p::Evaluate(good,sample.y,Vector::Constant(2,-1000)).valid);
    const auto zero=p::Fit(good,Vector::Zero(sample.y.size()),Vector::Constant(2,.5));
    EXPECT_NE(zero.at("runtime_convergence"),"passed");
    EXPECT_TRUE(zero.at("primary").at("kkt_passed").as_bool());
    EXPECT_EQ(zero.at("primary").at("rss"),0.0);
    EXPECT_EQ(zero.at("width_spectrum").at("rank"),0);
    // Geometry and observations are immutable across all evaluations.
    EXPECT_EQ(good.rows,sample.y.size()); EXPECT_EQ(good.atoms[0].size(),p::Domain(sample.grid,sample.atoms).atoms[0].size());
}

TEST(JointComponentNumericsTest, TrueWidthDoubleAndFloat32LinearControls)
{
    const Sample sample; const p::Domain domain(sample.grid,sample.atoms);
    const Vector eta=Eigen::Vector2d(.42,.67).array().log();
    const Vector truth=Eigen::Vector4d(2,-.3,1.5,.4);
    for(bool quantized:{false,true})
    {
        Vector y=sample.y;
        if(quantized) for(Eigen::Index k=0;k<y.size();++k) y(k)=static_cast<float>(y(k));
        const auto primary=p::Evaluate(domain,y,eta),reference=p::Evaluate(domain,y,eta,true);
        ASSERT_TRUE(primary.valid && reference.valid);
        EXPECT_TRUE(primary.certificate.at("kkt_passed").as_bool());
        EXPECT_LT((primary.beta-reference.beta).norm(),1e-10);
        EXPECT_LT((primary.beta-truth).norm(),quantized ? 1e-6 : 1e-12);
        if(quantized) EXPECT_GT(primary.residual.norm(),1e-9);
        else EXPECT_LT(primary.residual.norm(),1e-12);
    }
}

TEST(JointTestNumericsTest, TiledDerivativeMatchesDenseAcrossTileBoundaries)
{
    namespace n=p::runtime;
    Sample sample; const p::Domain base(sample.grid,sample.atoms);
    for(Eigen::Index rows:{8191,8192,8193})
    {
        std::vector<std::vector<p::Support>> support(base.atoms.size()); Vector y(rows);
        for(Eigen::Index first=0;first<rows;first+=base.rows)
            for(std::size_t a=0;a<support.size();++a) for(const auto & s:base.atoms[a])
                if(first+s.row<rows) support[a].push_back({first+s.row,s.square});
        for(Eigen::Index r=0;r<rows;++r) y(r)=sample.y(r%base.rows)+.03*std::sin(static_cast<double>(r));
        const p::Domain domain(rows,std::move(support)); const auto context=p::MakeContext(y,2);
        const auto e=p::Evaluate(domain,y,Eigen::Vector2d(.55,.51).array().log(),false,&context);
        ASSERT_TRUE(e.valid); const auto dense=p::DenseDifferentiate(e,context.scale,&context);
        ASSERT_TRUE(dense.valid); ASSERT_GT((dense.projected-dense.jacobian).norm(),1e-3);
        const auto correction=p::DenseLocalCorrection(e,dense,context);
        for(Eigen::Index tile:{17,8192})
        {
            n::DerivativeWorkForTesting()={};
            const auto d=n::PrepareDerivative(e,context.scale,&context,-1,tile);
            const auto reduced=n::ReduceDerivative(d,e.residual,true,tile);
            ASSERT_TRUE(reduced.valid);
            EXPECT_LE(n::DerivativeWorkForTesting().maximum_generated_rows,tile);
            EXPECT_LE(n::DerivativeWorkForTesting().maximum_reduction_rows,tile+e.x.cols());
            const auto actual=p::MaterializeDerivative(e,context.scale,&context,-1,tile);
            EXPECT_LT((actual.projected-dense.projected).norm()/dense.projected.norm(),1e-8);
            EXPECT_LT((actual.jacobian-dense.jacobian).norm()/dense.jacobian.norm(),1e-8);
            for(const auto & matrices:{std::pair{dense.projected,reduced.projected},std::pair{dense.jacobian,reduced.jacobian}})
            {
                const auto left=n::ComputeSpectrum(matrices.first,context.rank,2,false);
                const auto right=n::ComputeSpectrum(matrices.second,context.rank,2,false);
                EXPECT_EQ(left.rank,right.rank);
                EXPECT_LE((left.singular_values-right.singular_values).lpNorm<Eigen::Infinity>()/left.singular_values(0),1e-10);
                EXPECT_NEAR(left.threshold,right.threshold,1e-10*left.threshold);
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(reduced.jacobian,Eigen::ComputeThinU|Eigen::ComputeThinV);
            svd.setThreshold(context.rank.Relative(2));
            const Vector compact=svd.solve(-reduced.response);
            EXPECT_LE(((correction-compact).array().abs()/(1+correction.array().abs().max(compact.array().abs()))).maxCoeff(),1e-10);
            EXPECT_LT((reduced.projected_norms-dense.projected.colwise().norm().transpose()).norm()/dense.projected.norm(),1e-10);
            EXPECT_LT((reduced.jacobian_norms-dense.jacobian.colwise().blueNorm().transpose()).norm()/dense.jacobian.norm(),1e-10);
        }
    }
}

TEST(JointComponentNumericsTest, CompactJacobianGramIdentityDiagnosticsAcrossLattices)
{
    CompactJacobianSummary summary;
    for(const std::string topology:{"chain","cube"}) for(const int atoms:{8,32,128,256})
    {
        rhbm_gem::core::JointProblem problem(second_stage_test::OperatorWorkload(topology,atoms));
        const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
        const Vector eta=Vector::Constant(atoms,std::log(.55));
        const auto e=n::EvaluateProfile(data.domain,data.y,eta,false,&data.context);
        ASSERT_TRUE(e.valid)<<topology<<"-"<<atoms<<": "<<e.reason;
        CheckCompactJacobianIdentity(e,data.context,topology+"-"+std::to_string(atoms),true,true,true,summary);
    }
    const auto record=[](const char * key,double value) {
        std::ostringstream out; out<<std::setprecision(17)<<value; RecordProperty(key,out.str());
    };
    record("max_normalized_ZtP",summary.maximum.ztp_relative);
    record("max_normalized_Ztr",summary.maximum.ztr_relative);
    record("max_raw_normalized_Ztr",summary.maximum.ztr_relative_all);
    record("max_scaled_Ztr_stationarity",summary.maximum.ztr_scaled);
    record("max_full_gram_relative_frobenius_error",summary.maximum.gram_relative);
    record("max_full_gram_absolute_frobenius_error",summary.maximum.gram_absolute);
    record("max_full_gram_max_elementwise_error",summary.maximum.gram_elementwise);
    record("max_orthogonality_only_gram_relative_error",summary.maximum.orthogonal_gram_relative);
    record("max_compact_candidate_spectrum_relative_error",summary.maximum.spectrum_relative);
    record("max_weak_direction_projector_difference",summary.maximum.weak_projector_difference);
    record("max_response_coordinate_norm_difference",summary.maximum.response_norm_difference);
    record("max_response_gradient_difference",summary.maximum.response_gradient_difference);
    record("max_correction_difference",summary.maximum.correction_difference);
    RecordProperty("compact_candidate_rank_parity",summary.rank_parity ? "true" : "false");
}

TEST(JointComponentNumericsTest, CompactJacobianDiagnosticsCoverCancellationAndActiveFaces)
{
    CompactJacobianSummary cancellation;
    {
        Sample sample; const p::Domain domain(sample.grid,sample.atoms);
        const auto context=n::CreateContext(sample.y,2);
        auto e=n::EvaluateProfile(domain,sample.y,Eigen::Vector2d(.55,.51).array().log(),false,&context);
        ASSERT_TRUE(e.valid);
        e.derivative=e.x;
        n::SparseWorkForTesting()={};
        const auto prepared=n::PrepareDerivative(e,context.scale,&context);
        ASSERT_TRUE(prepared.valid);
        EXPECT_EQ(prepared.reference_order,n::SparseBackendEnabled());
        CheckCompactJacobianIdentity(e,context,"near-complete-cancellation",false,true,false,cancellation);
    }
    CompactJacobianSummary active_faces; bool saw_active=false,saw_interior=false;
    {
        Sample sample; sample.atoms.resize(1); sample.atoms[0].amplitude=.05; sample.atoms[0].charge=-1;
        for(Eigen::Index k=0;k<sample.y.size();++k)
            sample.y(k)=m::unique_grid::Direct(sample.grid.voxels[static_cast<std::size_t>(k)].position,sample.atoms,2.5);
        const p::Domain domain(sample.grid,sample.atoms); const auto context=n::CreateContext(sample.y,1);
        for(const double width:{.25,.42,.8})
        {
            const Vector eta=Vector::Constant(1,std::log(width));
            const auto e=n::EvaluateProfile(domain,sample.y,eta,false,&context);
            ASSERT_TRUE(e.valid)<<width<<": "<<e.reason;
            saw_active|=e.beta(0)==0; saw_interior|=e.beta(0)>0;
            CheckCompactJacobianIdentity(e,context,"active-face-width-"+std::to_string(width),true,true,false,active_faces);
        }
    }
    EXPECT_TRUE(saw_active && saw_interior);
    const auto record=[](const char * key,double value) {
        std::ostringstream out; out<<std::setprecision(17)<<value; RecordProperty(key,out.str());
    };
    record("cancellation_ZtP_absolute",cancellation.maximum.ztp_absolute);
    record("cancellation_full_gram_relative_error",cancellation.maximum.gram_relative);
    RecordProperty("cancellation_candidate_rank_parity",cancellation.rank_parity ? "true" : "false");
    record("active_face_max_normalized_ZtP",active_faces.maximum.ztp_relative);
    record("active_face_max_normalized_Ztr",active_faces.maximum.ztr_relative);
    record("active_face_max_raw_normalized_Ztr",active_faces.maximum.ztr_relative_all);
    record("active_face_max_scaled_Ztr_stationarity",active_faces.maximum.ztr_scaled);
    RecordProperty("active_face_candidate_rank_parity",active_faces.rank_parity ? "true" : "false");
    RecordProperty("active_and_interior_faces_seen",saw_active && saw_interior ? "true" : "false");
}

TEST(JointComponentNumericsTest, CompactJacobianRankBoundariesRespectPolicyAndFallBackWithinOneUlp)
{
    const n::RankPolicy policy{32,2,1};
    const double relative=policy.Relative(2);
    const std::array<n::RankBoundary,2> boundaries{n::RankBoundary::SvdNative,n::RankBoundary::StrictGreater};
    const std::array<double,3> scales{.125,1.,8.};
    for(const auto boundary:boundaries) for(const double scale:scales) for(const bool permute_rows:{false,true})
    {
        const double maximum=2.*scale,cutoff=relative*maximum;
        const std::array<std::pair<const char *,double>,5> grid{{
            {"clearly-below",.5*cutoff},
            {"just-below",std::nextafter(cutoff,0.)},
            {"boundary",cutoff},
            {"just-above",std::nextafter(cutoff,std::numeric_limits<double>::infinity())},
            {"clearly-above",2.*cutoff}}};
        for(const auto [label,minimum]:grid)
        {
            SCOPED_TRACE(std::string(label)+" scale="+std::to_string(scale)+
                " boundary="+std::to_string(static_cast<int>(boundary))+" permuted="+std::to_string(permute_rows));
            Matrix current=Matrix::Zero(2,2); current(0,0)=maximum; current(1,1)=minimum;
            Vector response(2); response<<.25*scale,-.125*scale;
            if(permute_rows) {current.row(0).swap(current.row(1)); std::swap(response(0),response(1));}
            Matrix compact=Matrix::Zero(4,2); compact.bottomRows(2)=current;
            Vector compact_response=Vector::Zero(4); compact_response.tail(2)=response;
            if(permute_rows)
            {
                Matrix reordered(4,2); Vector reordered_response(4);
                const std::array<Eigen::Index,4> order{3,0,2,1};
                for(Eigen::Index row=0;row<4;++row)
                {reordered.row(row)=compact.row(order[static_cast<std::size_t>(row)]); reordered_response(row)=compact_response(order[static_cast<std::size_t>(row)]);}
                compact=std::move(reordered); compact_response=std::move(reordered_response);
            }
            const n::RankRequest request{policy,2,-1.,boundary};
            const auto reference=n::EvaluateRank(current,request,&response,n::CompactSvdVectors::Right);
            const auto candidate=n::EvaluateRank(compact,request,&compact_response,n::CompactSvdVectors::Right);
            ASSERT_TRUE(reference.valid); ASSERT_TRUE(candidate.valid);
            ASSERT_EQ(reference.singular_values.size(),2); ASSERT_EQ(candidate.singular_values.size(),2);
            const std::string point=label;
            if(point=="clearly-below") EXPECT_EQ(reference.rank,1);
            if(point=="clearly-above") EXPECT_EQ(reference.rank,2);
            if(boundary==n::RankBoundary::StrictGreater)
                EXPECT_EQ(reference.rank,(reference.singular_values.array()>reference.threshold).count());
            else
            {
                Eigen::JacobiSVD<Matrix> native(current,Eigen::ComputeFullU|Eigen::ComputeFullV);
                native.setThreshold(relative); EXPECT_EQ(reference.rank,native.rank());
            }
            EXPECT_DOUBLE_EQ(reference.threshold,relative*maximum);
            EXPECT_DOUBLE_EQ(candidate.threshold,reference.threshold);
            EXPECT_LE((reference.singular_values-candidate.singular_values).norm()/maximum,1e-14);
            const bool boundary_ambiguous=reference.rank!=candidate.rank;
            if(!boundary_ambiguous)
            {
                const std::string reference_reason=reference.rank==2 ? "rank-verified-full" : "rank-deficient";
                const std::string candidate_reason=candidate.rank==2 ? "rank-verified-full" : "rank-deficient";
                EXPECT_EQ(reference_reason,candidate_reason);
            }
            const Matrix reference_projector=reference.right_vectors.col(1)*reference.right_vectors.col(1).transpose();
            const Matrix candidate_projector=candidate.right_vectors.col(1)*candidate.right_vectors.col(1).transpose();
            EXPECT_LE((reference_projector-candidate_projector).norm(),1e-12);
            EXPECT_LE((reference.solution-candidate.solution).norm(),1e-12*(1.+reference.solution.norm()));

            const double below=std::nextafter(reference.threshold,0.);
            const double above=std::nextafter(reference.threshold,std::numeric_limits<double>::infinity());
            const bool near_boundary=minimum>=below && minimum<=above;
            const bool expected_near=std::string(label)=="just-below" || std::string(label)=="boundary" ||
                std::string(label)=="just-above";
            EXPECT_EQ(near_boundary,expected_near);
            const bool compact_safe=n::CompactRankDecisionSafe(candidate);
            EXPECT_EQ(compact_safe,!expected_near);
            if(boundary_ambiguous) EXPECT_FALSE(compact_safe);
            const auto & selected=(compact_safe && !boundary_ambiguous) ? candidate : reference;
            EXPECT_EQ(selected.rank,reference.rank);
            EXPECT_DOUBLE_EQ(selected.threshold,reference.threshold);
            EXPECT_EQ(!near_boundary,std::string(label)=="clearly-below" || std::string(label)=="clearly-above");
        }
    }
}


TEST(JointComponentNumericsTest, SparseWorkspaceChecksPatternScopePolicyAndGeneration)
{
    namespace n=p::runtime;
    if(!n::SparseBackendEnabled()) GTEST_SKIP()<<"Optional SPQR backend";
    Matrix dense(6,2); dense<<1,0,2,1,0,3,4,2,1,1,0,2;
    n::Sparse x=dense.sparseView(); n::LinearWorkspace workspace;
    n::LinearPolicy policy{1e-14}; int domain{},other{};
    workspace.Bind(&domain,&domain,&policy); n::SparseWorkForTesting()={};
    auto first=workspace.Factor(x,{0,1},1e-14);
    ASSERT_EQ(first->Rank(),2);
    const Matrix rhs=Matrix::Ones(6,1);
    EXPECT_LT((first->LeastSquares(rhs)-dense.colPivHouseholderQr().solve(rhs)).norm(),1e-12);
    Matrix compact=first->Compact();
    EXPECT_LT((compact.transpose()*compact-dense.transpose()*dense).norm(),1e-12);
    x.valuePtr()[0]+=.1;
    auto second=workspace.Factor(x,{0,1},1e-14);
    EXPECT_EQ(n::SparseWorkForTesting().symbolic,1);
    EXPECT_EQ(n::SparseWorkForTesting().numeric,2);
    EXPECT_EQ(n::SparseWorkForTesting().symbolic_reuses,1);
    EXPECT_FALSE(first->Matches(x,{0,1}));
    EXPECT_THROW(first->LeastSquares(rhs),std::logic_error);
    EXPECT_TRUE(second->Matches(x,{0,1}));
    x.coeffRef(0,1)=.2; x.makeCompressed();
    workspace.Factor(x,{0,1},1e-14);
    EXPECT_EQ(n::SparseWorkForTesting().symbolic,2);
    workspace.Factor(x,{2,3},1e-14);
    EXPECT_EQ(n::SparseWorkForTesting().symbolic,3);
    policy.release_factor=64; workspace.Bind(&domain,&domain,&policy);
    workspace.Factor(x,{2,3},1e-14);
    EXPECT_EQ(n::SparseWorkForTesting().symbolic,4);
    workspace.Bind(&other,&domain,&policy); workspace.Factor(x,{2,3},1e-14);
    EXPECT_EQ(n::SparseWorkForTesting().symbolic,5);
    x.valuePtr()[0]=0;
    workspace.Factor(x,{2,3},1e-14);
    EXPECT_EQ(n::SparseWorkForTesting().symbolic,5); // Stored zero does not change CSC pattern.
    x.prune(0.); workspace.Factor(x,{2,3},1e-14);
    EXPECT_EQ(n::SparseWorkForTesting().symbolic,6);
    workspace.Factor(x,{2,3},0.);
    EXPECT_EQ(n::SparseWorkForTesting().symbolic,7);
}

TEST(JointComponentNumericsTest, ExactNumericReuseCensusRequiresAllFactorInputsToMatch)
{
    namespace n=p::runtime;
    if(!n::SparseBackendEnabled()) GTEST_SKIP()<<"Optional SPQR backend";
    Matrix dense(6,2); dense<<1,0,2,1,0,3,4,2,1,1,0,2;
    n::Sparse x=dense.sparseView(); x.makeCompressed(); n::LinearWorkspace workspace;
    n::LinearPolicy policy{1e-14}; int domain{};
    n::SparseWorkForTesting()={}; n::LinearTelemetryScopeForTesting telemetry(true);
    n::ProfileEvaluationRoleScopeForTesting role(n::ProfileEvaluationRole::InitialProfile);
    workspace.Bind(&domain,&domain,&policy);
    workspace.Factor(x,{0,1},1e-14);
    workspace.Factor(x,{0,1},1e-14);
    x.valuePtr()[0]+=.1;
    workspace.Factor(x,{0,1},1e-14);
    workspace.Factor(x,{1,0},1e-14);
    auto changed=policy; changed.release_factor=64;
    workspace.Bind(&domain,&domain,&changed);
    workspace.Factor(x,{0,1},1e-14);

    const auto census=n::SparseWorkForTesting();
    EXPECT_EQ(census.numeric_factor_requests,5u);
    EXPECT_EQ(census.numeric_factor_exact_reuse_opportunities,1u);
    EXPECT_EQ(census.initial_profile_exact_reuse_opportunities,1u);
    EXPECT_EQ(census.trial_profile_exact_reuse_opportunities,0u);
    EXPECT_EQ(census.reference_exact_reuse_opportunities,0u);
    EXPECT_EQ(census.accepted_endpoint_exact_reuse_opportunities,0u);
    EXPECT_EQ(census.numeric_factor_pattern_only_matches,1u);
    EXPECT_EQ(census.numeric_factor_value_mismatches,1u);
    EXPECT_EQ(census.numeric_factor_column_mismatches,1u);
    EXPECT_EQ(census.numeric_factor_policy_mismatches,1u);
    EXPECT_EQ(census.numeric_factor_pattern_mismatches,0u);

    n::SparseWorkForTesting()={}; n::LinearWorkspace fresh;
    fresh.Bind(&domain,&domain,&policy); fresh.Factor(x,{0,1},1e-14);
    EXPECT_EQ(n::SparseWorkForTesting().numeric_factor_requests,1u);
    EXPECT_EQ(n::SparseWorkForTesting().numeric_factor_exact_reuse_opportunities,0u);
}

TEST(JointComponentNumericsTest, SparseWeightedSolveAndIndependentReferenceMatchDenseOracle)
{
    namespace n=p::runtime;
    Matrix x(9,4); x<<1,2,0,0,2,0,1,0,0,3,0,1,1,0,4,0,0,1,2,3,2,0,0,1,1,1,1,1,3,2,1,0,0,0,0,0;
    const Vector truth=Eigen::Vector4d(2,-.3,1.5,.4);
    Vector y=x*truth; y(4)+=.03; y(8)=2;
    Vector w(9); w<<1,.5,0,2,1,3,.2,1,1;
    const n::Sparse sparse=x.sparseView();
    const auto primary=n::SolveLinear(sparse,y,w),reference=n::SolveLinear(sparse,y,w,true),dense=n::SolveLinear(x,y,w,true);
    ASSERT_TRUE(primary.valid && reference.valid && dense.valid);
    EXPECT_EQ(primary.rank,dense.rank); EXPECT_EQ(reference.rank,dense.rank);
    EXPECT_LT((primary.beta-dense.beta).norm(),1e-10);
    EXPECT_LT((reference.beta-dense.beta).norm(),1e-10);
    if(n::SparseBackendEnabled())
    {
        const auto before=n::SparseWorkForTesting().reference;
        const auto numeric=n::SparseWorkForTesting().numeric;
        ASSERT_TRUE(primary.factor);
        n::SolveLinear(sparse,y,w,true);
        EXPECT_EQ(n::SparseWorkForTesting().reference,before+1);
        EXPECT_EQ(n::SparseWorkForTesting().numeric,numeric);
        EXPECT_EQ(primary.factor->Rank(),primary.rank);
    }
}

TEST(JointComponentNumericsTest, SparseDerivativeReusesOnlyTheCanonicalCurrentFace)
{
    namespace n=p::runtime;
    if(!n::SparseBackendEnabled()) GTEST_SKIP()<<"Optional SPQR backend";
    Sample sample; const p::Domain domain(sample.grid,sample.atoms);
    const auto context=p::MakeContext(sample.y,2); const Vector eta=Eigen::Vector2d(.55,.51).array().log();
    n::LinearWorkspace workspace;
    auto e=n::EvaluateProfile(domain,sample.y,eta,false,&context,nullptr,&workspace);
    ASSERT_TRUE(e.valid); ASSERT_TRUE(e.factor);
    n::SparseWorkForTesting()={};
    ASSERT_TRUE(n::PrepareDerivative(e,context.scale,&context).valid);
    EXPECT_EQ(n::SparseWorkForTesting().factor_reuses,1);
    EXPECT_EQ(n::SparseWorkForTesting().numeric,0);
    e.beta(0)=0; // A free-set flag is not the canonical derivative face.
    ASSERT_TRUE(n::PrepareDerivative(e,context.scale,&context).valid);
    EXPECT_EQ(n::SparseWorkForTesting().numeric,1);
    const auto raw=n::EvaluateState(domain,sample.y,eta,e.beta,context);
    const auto copy=raw.beta;
    ASSERT_TRUE(n::PrepareDerivative(raw,context.scale,&context).valid);
    EXPECT_EQ(raw.beta,copy);
}

TEST(JointComponentNumericsTest, SparseReferencePreservesSvdRankNearDegeneracy)
{
    namespace n=p::runtime;
    for(double delta:{1e-8,1e-12,1e-15,0.})
    {
        Matrix x=Matrix::Zero(9,4);
        x(0,0)=1; x(1,1)=1; x(2,2)=1; x(0,3)=1; x(3,3)=delta;
        const Vector y=x*Eigen::Vector4d(2,.2,1,.3),w=Vector::Ones(9);
        const n::Sparse sparse=x.sparseView();
        const auto reference=n::SolveLinear(sparse,y,w,true),dense=n::SolveLinear(x,y,w,true);
        EXPECT_EQ(reference.rank,dense.rank)<<delta;
        EXPECT_EQ(reference.valid,dense.valid)<<delta;
        if(reference.valid) EXPECT_LT((x*reference.beta-x*dense.beta).norm(),1e-10);
    }
}

TEST(JointComponentNumericsTest, SparseCancellationRetainsTiledDerivativePrecision)
{
    namespace n=p::runtime;
    if(!n::SparseBackendEnabled()) GTEST_SKIP()<<"Optional SPQR backend";
    Sample sample; const p::Domain domain(sample.grid,sample.atoms);
    const auto context=p::MakeContext(sample.y,2);
    auto e=n::EvaluateProfile(domain,sample.y,Eigen::Vector2d(.55,.51).array().log(),false,&context);
    ASSERT_TRUE(e.valid);
    e.derivative=e.x; // Each raw width column now lies in the free design span.
    n::SparseWorkForTesting()={};
    const auto prepared=n::PrepareDerivative(e,context.scale,&context);
    ASSERT_TRUE(prepared.valid);
    EXPECT_EQ(n::SparseWorkForTesting().cancellation_reductions,1);
    const auto dense=p::DenseDifferentiate(e,context.scale,&context);
    Matrix projected,jacobian; prepared.Rows(0,e.x.rows(),projected,jacobian);
    EXPECT_LT((projected-dense.projected).norm(),1e-12);
    EXPECT_LT((jacobian-dense.jacobian).norm(),1e-12);
}

TEST(JointComponentNumericsTest, CompactSvdDispatchAndRepeatedSpectrumMatchJacobi)
{
    namespace n=p::runtime;
    for(Eigen::Index size:{15,16,17,32})
    {
        Matrix x(size+3,size);
        for(Eigen::Index i=0;i<x.rows();++i) for(Eigen::Index k=0;k<size;++k)
            x(i,k)=std::sin(static_cast<double>(i*size+k+1));
        const Eigen::HouseholderQR<Matrix> qr(x);
        const Matrix q=qr.householderQ()*Matrix::Identity(x.rows(),size);
        Vector values=Vector::Ones(size); values.tail(size/2).setConstant(.25);
        x=q*values.asDiagonal();
        const Vector rhs=Vector::LinSpaced(x.rows(),-.3,2.);
        const double relative=std::numeric_limits<double>::epsilon()*1000000;
        n::CompactSvdResult expected;
        {SvdMode mode(n::CompactSvdMode::Legacy); expected=n::CompactSvd(x,relative,-1,&rhs);}
        n::SparseWorkForTesting()={};
        const auto spectrum=n::CompactSvd(x,relative),solved=n::CompactSvd(x,relative,-1,&rhs);
        ASSERT_TRUE(expected.valid && spectrum.valid && solved.valid);
        EXPECT_EQ(spectrum.used_bdc,size>=16); EXPECT_FALSE(spectrum.jacobi_retry);
        EXPECT_EQ(spectrum.rank,size); EXPECT_EQ(solved.rank,expected.rank);
        EXPECT_LT((spectrum.singular_values-expected.singular_values).lpNorm<Eigen::Infinity>(),1e-10);
        EXPECT_LT((solved.solution-expected.solution).norm(),1e-10);
        EXPECT_NEAR(spectrum.threshold,relative,1e-10*relative);
        EXPECT_EQ(n::SparseWorkForTesting().bdc_svds,size>=16 ? 2 : 0);
        EXPECT_EQ(n::SparseWorkForTesting().reference_solves,1);
        EXPECT_EQ(n::SparseWorkForTesting().free_design_svds,1);
    }
}

TEST(JointComponentNumericsTest, CompactSvdPreservesRankBoundaryAndAbsoluteOverride)
{
    namespace n=p::runtime;
    constexpr double relative=1e-9;
    for(double absolute:{-1.,1e-7}) for(double factor:{0.,.5,.999999,1.,1.000001,2.})
    {
        const double threshold=absolute<0 ? relative : absolute;
        Matrix x=Matrix::Identity(16,16); x(15,15)=threshold*factor;
        const Vector rhs=x*Vector::LinSpaced(16,.5,2.);
        n::CompactSvdResult expected;
        {SvdMode mode(n::CompactSvdMode::Legacy); expected=n::CompactSvd(x,relative,absolute,&rhs);}
        n::SparseWorkForTesting()={};
        const auto actual=n::CompactSvd(x,relative,absolute,&rhs);
        ASSERT_TRUE(actual.valid); EXPECT_EQ(actual.rank,expected.rank);
        EXPECT_EQ(actual.threshold,expected.threshold);
        EXPECT_LT((actual.solution-expected.solution).norm(),1e-10);
        const bool near=std::abs(factor-1)*threshold<=64*std::numeric_limits<double>::epsilon()*16;
        EXPECT_EQ(actual.jacobi_retry,near);
        EXPECT_EQ(n::SparseWorkForTesting().jacobi_retries,near ? 1 : 0);
        EXPECT_EQ(n::SparseWorkForTesting().reference_solves,1);
        if(near) EXPECT_GT(n::SparseWorkForTesting().jacobi_retry_seconds,0);
    }
}

TEST(JointComponentNumericsTest, CompactSvdRightVectorsPreserveCovarianceAndBoundaryRetry)
{
    namespace n=p::runtime;
    for(Eigen::Index size:{3,17,32})
    {
        Matrix x=Matrix::Zero(size+2,size);
        for(Eigen::Index k=0;k<size;++k)
        {x(k,k)=1+.03*static_cast<double>(k); x(size,k)=.02*static_cast<double>(k);}
        const Eigen::JacobiSVD<Matrix> reference(x,Eigen::ComputeThinV);
        const auto right=n::CompactSvd(x,1e-12,-1,nullptr,n::CompactSvdVectors::Right);
        ASSERT_TRUE(right.valid); ASSERT_EQ(right.right_vectors.cols(),size);
        EXPECT_EQ(right.solution.size(),0);
        EXPECT_EQ(n::CompactSvd(x,1e-12).right_vectors.size(),0);
        const Matrix factor=right.right_vectors*right.singular_values.cwiseInverse().asDiagonal();
        const Matrix expected=reference.matrixV()*reference.singularValues().cwiseInverse().asDiagonal();
        EXPECT_LT((factor*factor.transpose()-expected*expected.transpose()).norm(),1e-10);
        x.setZero(); x.topRows(size).setIdentity(); x(size-1,size-1)=1e-9;
        const auto boundary=n::CompactSvd(x,1e-9,-1,nullptr,n::CompactSvdVectors::Right);
        EXPECT_TRUE(boundary.valid); EXPECT_EQ(boundary.jacobi_retry,size>=16);
        EXPECT_TRUE(boundary.right_vectors.allFinite());
    }
}

TEST(JointComponentNumericsTest, InPlaceTiledQrMatchesOriginalOrthogonalTransforms)
{
    namespace n=p::runtime;
    for(Eigen::Index responses:{0,2})
    {
        n::TiledQR actual(4,responses); Matrix r(0,4),target(0,responses);
        for(Eigen::Index count:{3,5,3})
        {
            Matrix rows(count,4),rhs(count,responses);
            for(Eigen::Index i=0;i<rows.size();++i) rows.data()[i]=std::sin(static_cast<double>(i+count));
            for(Eigen::Index i=0;i<rhs.size();++i) rhs.data()[i]=std::cos(static_cast<double>(i+count));
            // Frozen copying implementation: compare R and transformed RHS,
            // including an underdetermined first tile and zero RHS columns.
            Matrix a(r.rows()+count,4),b(target.rows()+count,responses);
            a.topRows(r.rows())=r; a.bottomRows(count)=rows;
            b.topRows(target.rows())=target; b.bottomRows(count)=rhs;
            const Eigen::HouseholderQR<Matrix> qr(a);
            const Matrix transformed=qr.householderQ().adjoint()*b;
            const auto keep=std::min(a.rows(),a.cols());
            r=qr.matrixQR().topRows(keep).triangularView<Eigen::Upper>(); target=transformed.topRows(keep);
            actual.Append(rows,rhs);
            EXPECT_EQ((actual.r-r).norm(),0); EXPECT_EQ((actual.target-target).norm(),0);
        }
    }
}

TEST(JointComponentNumericsTest, CompactSvdZeroAndInvalidInputsCannotQualify)
{
    namespace n=p::runtime;
    Matrix x=Matrix::Zero(16,16); const Vector rhs=Vector::Ones(16);
    const auto zero=n::CompactSvd(x,1e-12,-1,&rhs);
    ASSERT_TRUE(zero.valid); EXPECT_EQ(zero.rank,0); EXPECT_TRUE(zero.jacobi_retry);
    EXPECT_EQ(zero.solution.norm(),0);
    x(0,0)=std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(n::CompactSvd(x,1e-12).valid);
    EXPECT_FALSE(n::CompactSvd(Matrix(0,0),1e-12).valid);
    EXPECT_FALSE(n::CompactSvd(Matrix::Identity(17,17),1e-12,-1,&rhs).valid);
}

TEST(JointComponentNumericsTest, CompactSvdLargeDerivativeAndReferencePreserveResidualCorrection)
{
    namespace n=p::runtime;
    Sample sample; const p::Domain base(sample.grid,sample.atoms);
    constexpr Eigen::Index copies=4;
    std::vector<std::vector<p::Support>> support(static_cast<std::size_t>(2*copies));
    Vector y(base.rows*copies),eta(2*copies);
    for(Eigen::Index k=0;k<copies;++k)
    {
        y.segment(k*base.rows,base.rows)=sample.y;
        eta.segment(2*k,2)=Eigen::Vector2d(.55,.51).array().log();
        for(std::size_t a=0;a<base.atoms.size();++a) for(const auto & s:base.atoms[a])
            support[static_cast<std::size_t>(2*k)+a].push_back({k*base.rows+s.row,s.square});
    }
    const p::Domain domain(y.size(),std::move(support)); const auto context=p::MakeContext(y,eta.size());
    const auto e=n::EvaluateProfile(domain,y,eta,false,&context);
    ASSERT_TRUE(e.valid); ASSERT_GT(e.residual.norm(),1e-3);
    n::TiledDifferential legacy; n::Evaluation reference;
    {SvdMode mode(n::CompactSvdMode::Legacy);
        legacy=n::PrepareDerivative(e,context.scale,&context);
        reference=n::EvaluateProfile(domain,y,eta,true,&context);}
    n::SparseWorkForTesting()={};
    const auto actual=n::PrepareDerivative(e,context.scale,&context);
    const auto fast=n::EvaluateProfile(domain,y,eta,true,&context);
    ASSERT_TRUE(actual.valid && legacy.valid && reference.valid && fast.valid);
    EXPECT_GT(n::SparseWorkForTesting().bdc_svds,0);
    EXPECT_LT((actual.coefficients-legacy.coefficients).norm(),1e-10);
    EXPECT_LT((actual.correction-legacy.correction).norm(),1e-10);
    EXPECT_LT((fast.beta-reference.beta).norm(),1e-10);
    EXPECT_TRUE(n::CheckTrust(domain,y,e,context,fast).passed);
    const auto dense=p::DenseDifferentiate(p::Evaluation(e),context.scale,&context);
    Matrix projected,jacobian; actual.Rows(0,y.size(),projected,jacobian);
    EXPECT_GT((projected-jacobian).norm(),1e-3);
    EXPECT_LT((jacobian-dense.jacobian).norm()/dense.jacobian.norm(),1e-8);
    EXPECT_LT((projected-dense.projected).norm()/dense.projected.norm(),1e-8);
    const auto rejected=n::PrepareDerivative(e,context.scale,&context,100.);
    EXPECT_FALSE(rejected.valid); EXPECT_EQ(rejected.reason,"rank-deficient-free-design");
    EXPECT_EQ(n::SparseWorkForTesting().derivative_preparations,2);
    EXPECT_GT(n::SparseWorkForTesting().derivative_seconds,0);
}

TEST(JointComponentNumericsTest, CompactSvdReferencePreservesWeightedChangingActiveFace)
{
    namespace n=p::runtime;
    Matrix x=Matrix::Zero(33,16); x.topRows(16).setIdentity(); x.row(16).setConstant(.2);
    Vector truth=Vector::Ones(16); truth(0)=-.5; truth(3)=-.3;
    const Vector y=x*truth; Vector weights=Vector::LinSpaced(33,.5,2.); weights(32)=0;
    const n::Sparse sparse=x.sparseView();
    const auto oracle=n::SolveLinear(x,y,weights,true);
    n::SparseWorkForTesting()={};
    const auto candidate=n::SolveLinear(sparse,y,weights,true,true);
    ASSERT_TRUE(oracle.valid && candidate.valid);
    EXPECT_EQ(candidate.beta(0),0); EXPECT_GT(n::SparseWorkForTesting().bdc_svds,0);
    EXPECT_GT(n::SparseWorkForTesting().reference_svds,1);
    EXPECT_EQ(candidate.rank,oracle.rank);
    EXPECT_LT((candidate.beta-oracle.beta).norm(),1e-10);
}

TEST(JointComponentNumericsTest, SharedRankContractPreservesBothBoundaryRules)
{
    namespace n=p::runtime;
    const double relative=n::RankPolicy{1000,0,0}.Relative(4);
    for(double absolute:{-1.,1e-7}) for(double multiple:{0.,.5,1.,2.})
    {
        const double threshold=absolute<0 ? relative : absolute;
        Matrix x=Matrix::Identity(4,4); x(3,3)=threshold*multiple;
        const auto old=n::CompactSvd(x,relative,absolute);
        const auto native=n::EvaluateRank(x,{{1000,8,4},4,absolute});
        const auto strict=n::EvaluateRank(x,{{1000,8,4},4,absolute,n::RankBoundary::StrictGreater});
        ASSERT_TRUE(native.valid && strict.valid);
        EXPECT_EQ(native.rank,old.rank); EXPECT_EQ(native.threshold,old.threshold);
        EXPECT_EQ(strict.rank,(old.singular_values.array()>old.threshold).count());
    }
    const auto global=n::EvaluateRank(Matrix::Identity(2,2),{{1000,4,2},4});
    ASSERT_TRUE(global.valid); EXPECT_EQ(global.threshold,n::RankPolicy({1000,4,2}).Relative(4));
    EXPECT_EQ(n::EvaluateRank(Matrix::Zero(4,4),{{1000,8,4},4}).rank,0);
    EXPECT_FALSE(n::EvaluateRank(Matrix::Constant(4,4,n::unavailable),{{1000,8,4},4}).valid);
}
