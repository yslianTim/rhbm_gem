#include <gtest/gtest.h>
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/joint_component/ResourceWork.hpp"
#include "core/detail/joint_component/LinearSolve.hpp"
#include "core/detail/joint_component/CompactSvd.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "support/JointTestNumerics.hpp"
#include <unsupported/Eigen/NonLinearOptimization>
#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <numbers>
#include <sstream>
#include <tuple>

namespace {
namespace m=second_stage_test::matched;
namespace p=m::joint_abc;
namespace n=p::runtime;
using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;
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

TEST(JointComponentNumericsTest, TiledDerivativePreservesCancellationReferenceFallback)
{
    std::vector<std::vector<n::Support>> support(1);
    support[0]={{0,0.0},{1,0.25},{2,0.25}};
    Vector y(3);
    for(Eigen::Index row=0;row<3;++row)
    {
        const double square=row==0 ? 0.0 : 0.25;
        const auto basis=n::EvaluateKernel(square,.55,2.5);
        y(row)=3.0*basis.gaussian+.25*basis.charge+(row==1 ? 1e-4 : row==2 ? -1e-4 : 0.0);
    }
    const n::Domain domain(3,std::move(support));
    const auto context=p::MakeContext(y,1);
    const Vector eta=Vector::Constant(1,std::log(.55));
    const auto evaluation=n::EvaluateProfile(domain,y,eta,&context);
    ASSERT_TRUE(evaluation.valid)<<evaluation.reason;
    ASSERT_GT(evaluation.beta(0),0.0);
    EXPECT_GT(evaluation.residual.norm(),0.0);

    const auto full=n::PrepareDerivative(evaluation,context.scale,&context,-1,3);
    ASSERT_TRUE(full.valid)<<full.reason;
    EXPECT_TRUE(full.reference_order);
    const auto dense=p::DenseDifferentiate(evaluation,context.scale,&context);
    ASSERT_TRUE(dense.valid)<<dense.reason;
    Matrix raw_tile,projected_tile,jacobian_tile;
    full.Rows(0,3,raw_tile,projected_tile,jacobian_tile);
    EXPECT_LT((projected_tile-dense.projected).norm()/dense.projected.norm(),1e-8);
    EXPECT_LT((jacobian_tile-dense.jacobian).norm()/dense.jacobian.norm(),1e-8);
    full.Rows(0,2,raw_tile,projected_tile,jacobian_tile);
    const auto * raw_storage=raw_tile.data();
    const auto * projected_storage=projected_tile.data();
    const auto * jacobian_storage=jacobian_tile.data();
    full.Rows(1,2,raw_tile,projected_tile,jacobian_tile);
    EXPECT_EQ(raw_tile.data(),raw_storage);
    EXPECT_EQ(projected_tile.data(),projected_storage);
    EXPECT_EQ(jacobian_tile.data(),jacobian_storage);
    const auto reduced=n::ReduceDerivative(full,evaluation.residual,3);
    ASSERT_TRUE(reduced.valid)<<reduced.reason;

    const auto projected=n::PrepareProjectedWidthDerivative(evaluation,context.scale,&context,-1,3);
    ASSERT_TRUE(projected.valid)<<projected.reason;
    EXPECT_TRUE(projected.reference_order);
    EXPECT_TRUE(projected.projected_only);
    Matrix raw,projected_only_tile,jacobian;
    projected.Rows(0,3,raw,projected_only_tile,jacobian);
    EXPECT_LT((projected_only_tile-dense.projected).norm()/dense.projected.norm(),1e-8);
    const auto projected_reduced=n::ReduceProjectedWidthDerivative(projected,3);
    ASSERT_TRUE(projected_reduced.valid)<<projected_reduced.reason;
}

TEST(JointComponentNumericsTest, BasisAssemblyMatchesTripletsForOrderedAndDuplicateSupports)
{
    auto expected=[](const n::Domain & domain,const Vector & eta) {
        n::Sparse x(domain.rows,2*eta.size()),derivative(domain.rows,2*eta.size());
        std::vector<Eigen::Triplet<double>> x_entries,derivative_entries;
        for(Eigen::Index a=0;a<eta.size();++a) for(const auto & support:domain.atoms[static_cast<std::size_t>(a)])
        {
            const auto basis=n::EvaluateKernel(support.square,std::exp(eta(a)),2.5);
            const std::array<double,4> values{basis.gaussian,basis.charge,
                basis.gaussian_log_width,basis.charge_log_width};
            for(Eigen::Index k=0;k<2;++k)
            {
                if(values[static_cast<std::size_t>(k)]!=0.0)
                    x_entries.emplace_back(support.row,2*a+k,values[static_cast<std::size_t>(k)]);
                if(values[static_cast<std::size_t>(k+2)]!=0.0)
                    derivative_entries.emplace_back(support.row,2*a+k,values[static_cast<std::size_t>(k+2)]);
            }
        }
        x.setFromTriplets(x_entries.begin(),x_entries.end());
        derivative.setFromTriplets(derivative_entries.begin(),derivative_entries.end());
        return std::pair{std::move(x),std::move(derivative)};
    };
    const Vector y=Vector::Ones(5),eta=Vector::Constant(1,std::log(.05));
    const std::vector<std::vector<n::Support>> ordered{{{0,0.0},{1,.0004},{2,.04},{3,6.25},{4,9.0}}};
    const std::vector<std::vector<n::Support>> unsorted_duplicates{{{2,.04},{0,0.0},{2,.01},
        {1,.0004},{3,6.25},{4,9.0}}};
    const std::vector<std::vector<n::Support>> noncontiguous{{{0,0.0},{2,.04},{4,9.0}}};
    const std::vector<std::vector<n::Support>> empty_support(1);
    for(const auto & support:{ordered,unsorted_duplicates,noncontiguous,empty_support})
    {
        const n::Domain domain(5,support);
        const auto context=p::MakeContext(y,1);
        const auto evaluation=n::EvaluateProfile(domain,y,eta,&context);
        const auto [expected_x,expected_derivative]=expected(domain,eta);
        EXPECT_EQ((evaluation.x-expected_x).norm(),0.0);
        EXPECT_EQ((evaluation.derivative-expected_derivative).norm(),0.0);
        EXPECT_EQ(evaluation.x.coeff(3,0),0.0);
        EXPECT_EQ(evaluation.derivative.coeff(3,0),0.0);
        EXPECT_EQ(evaluation.x.coeff(4,1),0.0);
        EXPECT_EQ(evaluation.derivative.coeff(4,1),0.0);
    }

    const n::Domain parent(5,{{{0,0.0},{2,.04},{3,.25}}});
    const auto local=parent.Select(n::Indices{0},3,
        n::CompactRowMapping::FromParentRows(n::Indices{3,0,2}));
    ASSERT_EQ(local.atoms[0][0].row,1);
    ASSERT_EQ(local.atoms[0][1].row,2);
    ASSERT_EQ(local.atoms[0][2].row,0);
    const Vector local_y=Vector::Ones(3);
    const auto context=p::MakeContext(local_y,1);
    const auto evaluation=n::EvaluateProfile(local,local_y,eta,&context);
    const auto [expected_x,expected_derivative]=expected(local,eta);
    EXPECT_EQ((evaluation.x-expected_x).norm(),0.0);
    EXPECT_EQ((evaluation.derivative-expected_derivative).norm(),0.0);
}

TEST(JointComponentNumericsTest, TiledQRReusesScratchForRepeatedAssembledShape)
{
    n::TiledQR qr(2,1);
    Matrix rows(2,2); rows<<1.0,2.0,3.0,5.0;
    Matrix rhs(2,1); rhs<<.25,-.5;
    qr.Append(rows,rhs);
    qr.Append(rows,rhs);
    const auto * row_storage=qr.scratch_rows.data();
    const auto * rhs_storage=qr.scratch_rhs.data();
    qr.Append(rows,rhs);
    EXPECT_EQ(qr.scratch_rows.data(),row_storage);
    EXPECT_EQ(qr.scratch_rhs.data(),rhs_storage);
}

TEST(JointComponentNumericsTest, SparseSolveCountsZeroDesignRowsWithNonzeroRhsInObjective)
{
    namespace n=p::runtime;
    Matrix dense=Matrix::Zero(9,4);
    dense(0,0)=1.0;
    dense(1,1)=1.0;
    dense(2,2)=1.0;
    dense(3,3)=1.0;
    dense(4,0)=.5;
    dense(4,2)=-.25;
    dense(4,3)=.125;
    const n::Sparse sparse=dense.sparseView(0.0,0.0);
    const Vector truth=(Vector(4)<<.8,-.3,1.2,.5).finished();
    Vector rhs=dense*truth;
    rhs(4)+=.125;
    rhs(5)=.5;
    rhs(6)=-.25;
    rhs(7)=.75;

    ASSERT_LT(sparse.nonZeros(),dense.size());
    n::LinearSolveWorkspace workspace;
    const auto solved=n::SolveLinear(sparse,rhs,nullptr,0,&workspace);
    ASSERT_TRUE(solved.valid)<<solved.reason;
    EXPECT_EQ(solved.rank,4);
    EXPECT_GE(solved.beta(0),0.0);
    EXPECT_GE(solved.beta(2),0.0);

    const Eigen::ColPivHouseholderQR<Matrix> reference(dense);
    ASSERT_EQ(reference.rank(),dense.cols());
    const Vector reference_beta=reference.solve(rhs);
    EXPECT_LE((solved.beta-reference_beta).norm()/(1+reference_beta.norm()),2e-8);

    const auto certificate=n::CertifyLinear(sparse,rhs,solved.beta);
    ASSERT_TRUE(certificate.available);
    EXPECT_TRUE(certificate.feasible);
    EXPECT_TRUE(certificate.kkt_passed);
    EXPECT_LE(certificate.projected_kkt,1e-10);
    const Vector residual=dense*solved.beta-rhs;
    EXPECT_NEAR(certificate.objective,.5*residual.squaredNorm(),1e-13);
    EXPECT_NEAR(certificate.objective,
        .5*(rhs.segment(5,3).squaredNorm()+residual.head(5).squaredNorm()),1e-13);
    EXPECT_GT(rhs.segment(5,3).squaredNorm(),0.0);
}

TEST(JointComponentNumericsTest, ControlledQrMatchesDenseReferenceAcrossPartialFinalTile)
{
    namespace n=p::runtime;
    constexpr Eigen::Index rows=2053;
    constexpr double pi=std::numbers::pi;
    const std::array<std::pair<Eigen::Index,double>,2> fixtures{{{4,.25},{24,1.0}}};
    for(const auto [columns,support_ratio]:fixtures)
    {
        Matrix dense=Matrix::Zero(rows,columns);
        const Eigen::Index supported=std::max<Eigen::Index>(1,
            static_cast<Eigen::Index>(std::llround(support_ratio*static_cast<double>(columns))));
        for(Eigen::Index row=0;row<rows;++row) for(Eigen::Index col=0;col<columns;++col)
        {
            const bool present=support_ratio==1.0 ||
                ((col-row%columns+columns)%columns)<supported;
            if(!present) continue;
            const double phase=2*pi*static_cast<double>((col+1)*(row+1))/static_cast<double>(rows);
            dense(row,col)=std::sin(phase)+.37*std::cos(phase)+.001*static_cast<double>(row+col);
        }
        Vector truth(columns);
        for(Eigen::Index col=0;col<columns;++col)
            truth(col)=col%2 ? -.12-.003*static_cast<double>(col) : .4+.02*static_cast<double>(col);
        Vector rhs=dense*truth;
        for(Eigen::Index row=0;row<rows;++row) rhs(row)+=.001*std::sin(.13*static_cast<double>(row));
        const n::Sparse sparse=dense.sparseView(0.0,0.0);
        const auto solved=n::SolveLinear(sparse,rhs);
        ASSERT_TRUE(solved.valid)<<"columns="<<columns<<" reason="<<solved.reason;
        EXPECT_EQ(solved.rank,columns);
        EXPECT_EQ(rows%1024,5);

        const Eigen::ColPivHouseholderQR<Matrix> reference(dense);
        ASSERT_EQ(reference.rank(),columns);
        const Vector reference_beta=reference.solve(rhs);
        EXPECT_LE((solved.beta-reference_beta).norm()/(1+reference_beta.norm()),2e-8);

        const auto certificate=n::CertifyLinear(sparse,rhs,solved.beta);
        ASSERT_TRUE(certificate.available);
        EXPECT_TRUE(certificate.feasible);
        EXPECT_TRUE(certificate.kkt_passed);
        EXPECT_LE(certificate.projected_kkt,1e-10);
        const Vector residual=dense*solved.beta-rhs;
        EXPECT_NEAR(certificate.objective,.5*residual.squaredNorm(),
            1e-12+2e-12*std::abs(certificate.objective));
    }
}

TEST(JointComponentNumericsTest, NearRankDeficientTwentyFourColumnDesignIsClassifiedAsRankTwentyThree)
{
    namespace n=p::runtime;
    constexpr Eigen::Index rows=2053,columns=24;
    constexpr double perturbation=1e-14;
    Matrix dense=Matrix::Zero(rows,columns);
    for(Eigen::Index col=0;col<columns-1;++col) dense(col,col)=1.0;
    dense(columns-3,columns-1)=1.0;
    dense(columns-2,columns-1)=perturbation;
    const n::Sparse sparse=dense.sparseView(0.0,0.0);
    Vector truth=Vector::LinSpaced(columns,.2,.8);
    const Vector rhs=dense*truth;
    const double rank_threshold=std::numeric_limits<double>::epsilon()*static_cast<double>(rows);
    Eigen::ColPivHouseholderQR<Matrix> reference(dense);
    reference.setThreshold(rank_threshold);
    ASSERT_EQ(reference.rank(),23);

    n::LinearPolicy policy;
    policy.rank_relative=rank_threshold;
    const auto solved=n::SolveLinear(sparse,rhs,&policy);
    EXPECT_FALSE(solved.valid);
    EXPECT_EQ(solved.reason,"rank-deficient");
    EXPECT_EQ(solved.rank,23);
    EXPECT_NE(solved.rank,columns);
}

;
;
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
        const auto ep=p::Evaluate(domain,sample.y,plus),em=p::Evaluate(domain,sample.y,minus);
        ASSERT_TRUE(ep.valid && em.valid);
        const Vector finite=(ep.residual-em.residual)/(2e-5*scale);
        EXPECT_LT((finite-d.jacobian.col(k)).norm()/finite.norm(),1e-7);
        EXPECT_NEAR(e.gradient(k),(ep.residual.squaredNorm()-em.residual.squaredNorm())/(4e-5*scale*scale),1e-9);
    }
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
        const auto plus=p::Evaluate(domain,sample.y,eta.array()+1e-5),minus=p::Evaluate(domain,sample.y,eta.array()-1e-5);
        ASSERT_EQ(plus.certificate.at("active_atoms"),minus.certificate.at("active_atoms"));
        EXPECT_LT(((plus.residual-minus.residual)/(2e-5*scale)-d.jacobian.col(0)).norm(),1e-8);
    }
    EXPECT_TRUE(active && interior);
}

TEST(JointTestNumericsTest, InvalidBStillFailsAndUninformativeWidthCanConverge)
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
    EXPECT_EQ(zero.at("runtime_convergence"),"passed");
    EXPECT_TRUE(zero.at("primary").at("kkt_passed").as_bool());
    EXPECT_EQ(zero.at("primary").at("rss"),0.0);
    EXPECT_TRUE(zero.at("numerical_convergence_passed").as_bool());
    // Geometry and observations are immutable across all evaluations.
    EXPECT_EQ(good.rows,sample.y.size()); EXPECT_EQ(good.atoms[0].size(),p::Domain(sample.grid,sample.atoms).atoms[0].size());
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
        const auto e=p::Evaluate(domain,y,Eigen::Vector2d(.55,.51).array().log(),&context);
        ASSERT_TRUE(e.valid); const auto dense=p::DenseDifferentiate(e,context.scale,&context);
        ASSERT_TRUE(dense.valid); ASSERT_GT((dense.projected-dense.jacobian).norm(),1e-3);
        const auto correction=p::DenseLocalCorrection(e,dense,context);
        for(Eigen::Index tile:{17,8192})
        {
            const auto d=n::PrepareDerivative(e,context.scale,&context,-1,tile);
            const auto reduced=n::ReduceDerivative(d,e.residual,tile);
            const auto width_derivative=n::PrepareProjectedWidthDerivative(e,context.scale,&context,-1,tile);
            const auto reduced_width=n::ReduceProjectedWidthDerivative(width_derivative,tile);
            ASSERT_TRUE(reduced.valid);
            ASSERT_TRUE(reduced_width.valid);
            const auto actual=p::MaterializeDerivative(e,context.scale,&context,-1,tile);
            EXPECT_LT((actual.projected-dense.projected).norm()/dense.projected.norm(),1e-8);
            EXPECT_LT((actual.jacobian-dense.jacobian).norm()/dense.jacobian.norm(),1e-8);
            for(const auto & matrices:{std::pair{dense.projected,reduced_width.projected},
                std::pair{dense.jacobian,reduced.jacobian}})
            {
                const auto left=p::DenseSpectrum(matrices.first,context.rank,2,false);
                const auto right=p::DenseSpectrum(matrices.second,context.rank,2,false);
                EXPECT_EQ(left.rank,right.rank);
                EXPECT_LE((left.singular_values-right.singular_values).lpNorm<Eigen::Infinity>()/left.singular_values(0),1e-10);
                EXPECT_NEAR(left.threshold,right.threshold,1e-10*left.threshold);
            }
            const Vector compact=p::JacobiLeastSquaresReference(
                reduced.jacobian,-reduced.response,context.rank);
            EXPECT_LE(((correction-compact).array().abs()/(1+correction.array().abs().max(compact.array().abs()))).maxCoeff(),1e-10);
            EXPECT_LT((reduced_width.projected_norms-dense.projected.colwise().norm().transpose()).norm()/dense.projected.norm(),1e-10);
            EXPECT_LT((reduced.jacobian_norms-dense.jacobian.colwise().blueNorm().transpose()).norm()/dense.jacobian.norm(),1e-10);
        }
    }
}

TEST(JointComponentNumericsTest, ProjectedWidthKeepsRankDeficiencyAsEvidenceWhileProfileDerivativeRejectsIt)
{
    namespace n=p::runtime;
    const std::vector<std::vector<n::Support>> supports{
        {{0,0.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}},
        {{0,0.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}}};
    const n::Domain domain(5,supports);
    const Vector y=Vector::Zero(5),eta=Vector::Constant(2,std::log(.7)),beta=Vector::Zero(4);
    const auto context=p::MakeContext(y,2);
    const auto state=n::EvaluateState(domain,y,eta,beta,context);
    ASSERT_TRUE(state.valid)<<state.reason;

    const auto projected=n::PrepareProjectedWidthDerivative(state,context.scale,&context);
    ASSERT_TRUE(projected.valid)<<projected.reason;
    ASSERT_TRUE(projected.free_design_rank);
    ASSERT_TRUE(projected.free_design_rank->available);
    const Matrix dense_design=Matrix(projected.free_design);
    Eigen::JacobiSVD<Matrix> reference(dense_design);
    reference.setThreshold(context.rank.Relative(4));
    EXPECT_EQ(projected.free_design_rank->rank,reference.rank());
    EXPECT_LT(reference.rank(),dense_design.cols());

    const auto full=n::PrepareDerivative(state,context.scale,&context);
    EXPECT_FALSE(full.valid);
    EXPECT_EQ(full.reason,"rank-deficient-free-design");
}

TEST(JointComponentNumericsTest, SparseConstrainedSolvePreservesInteriorBoundaryAndReleaseFaces)
{
    namespace n=p::runtime;
    Matrix x=Matrix::Zero(7,4);
    x(0,0)=1; x(1,1)=1; x(2,2)=1; x(3,3)=1;
    const n::Sparse sparse=x.sparseView();
    for(const Eigen::Vector4d truth:{Eigen::Vector4d(2,-.3,1.5,.4),Eigen::Vector4d(-.5,-.3,1.5,.4)})
    {
        const Vector y=x*truth; const auto solved=n::SolveLinear(sparse,y);
        ASSERT_TRUE(solved.valid);
        const auto certificate=n::CertifyLinear(sparse,y,solved.beta,std::max(1.0,y.norm()));
        ASSERT_TRUE(certificate.available); EXPECT_TRUE(certificate.feasible); EXPECT_TRUE(certificate.kkt_passed);
        EXPECT_EQ(solved.releases,0);
        if(truth(0)>0) {EXPECT_TRUE(certificate.active_atoms.empty()); EXPECT_NEAR(solved.beta(0),truth(0),1e-12);}
        else {ASSERT_EQ(certificate.active_atoms,std::vector<Eigen::Index>{0}); EXPECT_EQ(solved.beta(0),0);}
    }

    Matrix release_design(7,4);
    release_design<<0,0,-1,2, 0,1,-1,-1, -2,0,-2,-2, 2,-2,-1,0,
        0,0,2,2, 0,2,-1,0, -1,1,2,-1;
    const Vector release_y=(Vector(7)<<0,1,2,1,1,2,0).finished();
    const auto release_sparse=n::Sparse(release_design.sparseView());
    const auto released=n::SolveLinear(release_sparse,release_y);
    ASSERT_TRUE(released.valid); EXPECT_EQ(released.releases,1); EXPECT_GT(released.beta(0),0);
    EXPECT_EQ(released.beta(2),0);
    const auto released_certificate=n::CertifyLinear(release_sparse,release_y,released.beta,
        std::max(1.0,release_y.norm()));
    ASSERT_TRUE(released_certificate.available); EXPECT_TRUE(released_certificate.feasible);
    EXPECT_TRUE(released_certificate.kkt_passed);

    Matrix deficient=Matrix::Zero(7,4);
    deficient(0,0)=1; deficient(1,1)=1; deficient(0,2)=1; deficient(2,3)=1;
    const auto rejected=n::SolveLinear(n::Sparse(deficient.sparseView()),Vector::Ones(7));
    EXPECT_FALSE(rejected.valid); EXPECT_EQ(rejected.reason,"rank-deficient"); EXPECT_EQ(rejected.rank,3);
}

TEST(JointComponentNumericsTest, UnweightedLinearSolvePreservesRecordedNumericalBaseline)
{
    namespace n=p::runtime;
    const auto solve=[](const Matrix & x,const Vector & y) {
        const n::Sparse sparse=x.sparseView();
        auto result=n::SolveLinear(sparse,y);
        auto certificate=result.valid ? n::CertifyLinear(sparse,y,result.beta) : n::Certificate{};
        return std::pair{std::move(result),std::move(certificate)};
    };
    const auto relative_match=[](double actual,double expected,double tolerance) {
        return std::abs(actual-expected)<=tolerance*std::abs(expected);
    };

    Matrix identity=Matrix::Zero(7,4);
    identity(0,0)=1; identity(1,1)=1; identity(2,2)=1; identity(3,3)=1;
    Vector positive_truth(4); positive_truth<<2,-.3,1.5,.4;
    const auto [full_rank,full_rank_certificate]=solve(identity,identity*positive_truth);
    ASSERT_TRUE(full_rank.valid); EXPECT_EQ(full_rank.reason,"solved"); EXPECT_EQ(full_rank.rank,4);
    EXPECT_EQ(full_rank.active_set_iterations,1u); EXPECT_EQ(full_rank.releases,0);
    EXPECT_EQ(full_rank.sparse_structure_build_count,1u); EXPECT_EQ(full_rank.tile_mapping_build_count,1u);
    EXPECT_EQ(full_rank.free_face_changes,0u); EXPECT_EQ(full_rank.free_face_unchanged,0u);
    EXPECT_GT(full_rank.solver_total_seconds,0.0); EXPECT_GT(full_rank.sparse_qr_seconds,0.0);
    EXPECT_LE((full_rank.beta-positive_truth).norm(),1e-14);
    EXPECT_DOUBLE_EQ(full_rank_certificate.objective,0.0);
    EXPECT_EQ(full_rank.refinement_eligibility_count,1u); EXPECT_EQ(full_rank.refinement_trigger_count,0u);

    Vector boundary_truth(4); boundary_truth<<-.5,-.3,1.5,.4;
    const auto [boundary,boundary_certificate]=solve(identity,identity*boundary_truth);
    ASSERT_TRUE(boundary.valid); EXPECT_EQ(boundary.reason,"solved"); EXPECT_EQ(boundary.rank,3);
    EXPECT_EQ(boundary.active_set_iterations,2u); EXPECT_EQ(boundary.releases,0);
    EXPECT_EQ(boundary.sparse_structure_build_count,2u); EXPECT_EQ(boundary.tile_mapping_build_count,2u);
    EXPECT_EQ(boundary.free_face_changes,1u); EXPECT_EQ(boundary.free_face_unchanged,0u);
    EXPECT_DOUBLE_EQ(boundary.beta(0),0.0); EXPECT_NEAR(boundary.beta(1),-.3,1e-14);
    EXPECT_NEAR(boundary.beta(2),1.5,1e-14); EXPECT_NEAR(boundary.beta(3),.4,1e-14);
    EXPECT_NEAR(boundary_certificate.objective,.125,1e-14);
    EXPECT_NEAR(boundary_certificate.projected_kkt,2.7755575615628914e-17,1e-15);

    Matrix release_design(7,4);
    release_design<<0,0,-1,2, 0,1,-1,-1, -2,0,-2,-2, 2,-2,-1,0,
        0,0,2,2, 0,2,-1,0, -1,1,2,-1;
    const auto [released,release_certificate]=solve(release_design,(Vector(7)<<0,1,2,1,1,2,0).finished());
    ASSERT_TRUE(released.valid); EXPECT_EQ(released.reason,"solved"); EXPECT_EQ(released.rank,3);
    EXPECT_EQ(released.active_set_iterations,4u); EXPECT_EQ(released.releases,1);
    EXPECT_EQ(released.sparse_structure_build_count,4u); EXPECT_EQ(released.tile_mapping_build_count,4u);
    EXPECT_EQ(released.free_face_changes,3u); EXPECT_EQ(released.free_face_unchanged,0u);
    const Vector expected_release=(Vector(4)<<.038674033149171297,.28176795580110497,0.,-.18784530386740331).finished();
    EXPECT_LE((released.beta-expected_release).norm(),1e-14);
    EXPECT_NEAR(release_certificate.objective,4.834254143646409,1e-13);
    EXPECT_NEAR(release_certificate.projected_kkt,5.5511151231257827e-17,1e-15);

    Matrix deficient=Matrix::Zero(7,4);
    deficient(0,0)=1; deficient(1,1)=1; deficient(0,2)=1; deficient(2,3)=1;
    const auto [rejected,unavailable]=solve(deficient,Vector::Ones(7));
    EXPECT_FALSE(rejected.valid); EXPECT_EQ(rejected.reason,"rank-deficient"); EXPECT_EQ(rejected.rank,3);
    EXPECT_EQ(rejected.active_set_iterations,1u); EXPECT_EQ(rejected.solves,1);
    EXPECT_FALSE(unavailable.available);

    Matrix scaled=Matrix::Zero(7,4);
    scaled(0,0)=1e-6; scaled(1,1)=1e3; scaled(2,2)=1e-2; scaled(3,3)=1e4;
    Vector scaled_truth(4); scaled_truth<<2,-.3,1.5,.4;
    const auto [nonuniform,nonuniform_certificate]=solve(scaled,scaled*scaled_truth);
    ASSERT_TRUE(nonuniform.valid); EXPECT_EQ(nonuniform.rank,4); EXPECT_EQ(nonuniform.active_set_iterations,1u);
    EXPECT_LE((nonuniform.beta-scaled_truth).norm(),1e-14);
    EXPECT_DOUBLE_EQ(nonuniform_certificate.objective,0.0);

    Matrix large=Matrix::Zero(512,48); Vector large_truth(48);
    for(Eigen::Index column=0;column<large.cols();++column)
    {
        const Eigen::Index frequency=column/2+1;
        large_truth(column)=column%2 ? .15*std::sin(static_cast<double>(column)) : .2+.002*static_cast<double>(column);
        for(Eigen::Index row=0;row<large.rows();++row)
        {
            const double phase=2*std::numbers::pi*static_cast<double>(frequency*row)/static_cast<double>(large.rows());
            large(row,column)=column%2 ? std::cos(phase) : std::sin(phase);
        }
    }
    const auto [large_result,large_certificate]=solve(large,large*large_truth);
    ASSERT_TRUE(large_result.valid); EXPECT_EQ(large_result.rank,48);
    EXPECT_EQ(large_result.active_set_iterations,1u); EXPECT_EQ(large_result.releases,0);
    EXPECT_EQ(large_result.sparse_structure_build_count,1u); EXPECT_EQ(large_result.tile_mapping_build_count,1u);
    EXPECT_LE((large_result.beta-large_truth).norm(),1e-12);
    EXPECT_LT(large_certificate.objective,1e-27);
    EXPECT_LT(large_certificate.projected_kkt,2e-16);

    struct RefinementBaseline {double delta,a,c,b,d,objective,kkt; std::size_t successes,failures; int solves;};
    const std::array<RefinementBaseline,4> refinement_baselines{{
        {1e-8,1000000.001636615,-1000000.001636615,.39999999996198982,-.19999999997244269,
            2.9245477447186279e-20,7.0582428790544327e-11,1,0,5},
        {1e-10,999999.86949451349,-999999.86949451349,.40000000000915337,-.20000000001303836,
            1.2527628537315677e-20,1.9641177573248569e-11,1,0,2},
        {1e-12,999966.90945700917,-999966.90945700917,.40000000001083835,-.20000000004484161,
            2.4204223979307309e-20,8.9422913518433234e-11,1,0,2},
        {1e-14,1015378.31601107,-1015378.3160110699,.39999999977109302,-.19999999958936632,
            5.232877687860051e-19,9.3132257461547852e-10,0,1,1}}};
    for(const auto & baseline:refinement_baselines)
    {
        constexpr Eigen::Index rows=32;
        Matrix x(rows,4);
        for(Eigen::Index row=0;row<rows;++row)
        {
            const double k=static_cast<double>(row);
            const double a=std::sin(.23*k)+.2*std::cos(.09*k);
            const double b=std::cos(.17*k)+.3*std::sin(.31*k);
            const double c=std::sin(.37*k)+.1*std::cos(.41*k);
            const double d=std::cos(.29*k)-.2*std::sin(.13*k);
            x(row,0)=a; x(row,1)=a+baseline.delta*b; x(row,2)=c; x(row,3)=d;
        }
        Vector truth(4); truth<<1e6,-1e6,.4,-.2;
        const auto [solved,certificate]=solve(x,x*truth);
        ASSERT_TRUE(solved.valid); EXPECT_EQ(solved.reason,"solved"); EXPECT_EQ(solved.rank,4);
        EXPECT_EQ(solved.active_set_iterations,1u); EXPECT_EQ(solved.releases,0);
        EXPECT_EQ(solved.refinement_eligibility_count,1u); EXPECT_EQ(solved.refinement_trigger_count,1u);
        EXPECT_EQ(solved.refinement_success_count,baseline.successes);
        EXPECT_EQ(solved.refinement_failure_count,baseline.failures);
        EXPECT_EQ(solved.refinement_factorizations,1u); EXPECT_EQ(solved.solves,baseline.solves);
        EXPECT_TRUE(relative_match(solved.beta(0),baseline.a,1e-12));
        EXPECT_TRUE(relative_match(solved.beta(1),baseline.c,1e-12));
        EXPECT_TRUE(relative_match(solved.beta(2),baseline.b,1e-12));
        EXPECT_TRUE(relative_match(solved.beta(3),baseline.d,1e-12));
        EXPECT_TRUE(relative_match(certificate.objective,baseline.objective,1e-11));
        EXPECT_TRUE(relative_match(certificate.projected_kkt,baseline.kkt,1e-11));
        EXPECT_EQ(certificate.kkt_passed,baseline.failures==0);
    }
}

TEST(JointComponentNumericsTest, ActiveSetWorkspaceResetsAcrossFacesShapesAndFailures)
{
    namespace n=p::runtime;
    n::LinearSolveWorkspace workspace;
    Matrix small=Matrix::Zero(7,4);
    small(0,0)=1; small(1,1)=1; small(2,2)=1; small(3,3)=1;
    Vector truth(4); truth<<2,-.3,1.5,.4;
    const auto small_sparse=n::Sparse(small.sparseView());
    const Vector small_y=small*truth;
    const auto compare_bitwise=[](const Vector & left,const Vector & right) {
        if(left.size()!=right.size()) {ADD_FAILURE()<<"Vector sizes differ."; return;}
        EXPECT_EQ(std::memcmp(left.data(),right.data(),static_cast<std::size_t>(left.size())*sizeof(double)),0);
    };

    const auto fresh=n::SolveLinear(small_sparse,small_y);
    const auto reused=n::SolveLinear(small_sparse,small_y,nullptr,0,&workspace);
    ASSERT_TRUE(fresh.valid && reused.valid);
    compare_bitwise(fresh.beta,reused.beta);

    Vector boundary_truth(4); boundary_truth<<-.5,-.3,1.5,.4;
    const Vector boundary_y=small*boundary_truth;
    const auto fresh_boundary=n::SolveLinear(small_sparse,boundary_y);
    const auto reused_boundary=n::SolveLinear(small_sparse,boundary_y,nullptr,0,&workspace);
    ASSERT_TRUE(fresh_boundary.valid && reused_boundary.valid);
    compare_bitwise(fresh_boundary.beta,reused_boundary.beta);
    EXPECT_EQ(reused_boundary.beta(0),0.0);

    Matrix deficient=Matrix::Zero(7,4);
    deficient(0,0)=1; deficient(1,1)=1; deficient(0,2)=1; deficient(2,3)=1;
    const auto failed=n::SolveLinear(n::Sparse(deficient.sparseView()),Vector::Ones(7),nullptr,0,&workspace);
    EXPECT_FALSE(failed.valid); EXPECT_EQ(failed.reason,"rank-deficient");

    Matrix wider=Matrix::Zero(64,8); Vector wider_truth(8);
    for(Eigen::Index column=0;column<wider.cols();++column)
    {
        wider(column,column)=1;
        wider_truth(column)=column%2 ? -.1*static_cast<double>(column) : .2+.1*static_cast<double>(column);
    }
    const auto wider_sparse=n::Sparse(wider.sparseView()); const Vector wider_y=wider*wider_truth;
    const auto fresh_wider=n::SolveLinear(wider_sparse,wider_y);
    const auto reused_wider=n::SolveLinear(wider_sparse,wider_y,nullptr,0,&workspace);
    ASSERT_TRUE(fresh_wider.valid && reused_wider.valid);
    compare_bitwise(fresh_wider.beta,reused_wider.beta);

    const auto retry=n::SolveLinear(small_sparse,small_y,nullptr,0,&workspace);
    ASSERT_TRUE(retry.valid);
    compare_bitwise(fresh.beta,retry.beta);
    const auto columns_capacity=workspace.columns.capacity(),entries_capacity=workspace.entries.capacity(),
        response_capacity=workspace.response.capacity(),local_column_capacity=workspace.local_column.capacity(),
        tile_columns_capacity=workspace.tile_columns.capacity();
    const auto repeated=n::SolveLinear(small_sparse,small_y,nullptr,0,&workspace);
    ASSERT_TRUE(repeated.valid);
    EXPECT_EQ(workspace.columns.capacity(),columns_capacity);
    EXPECT_EQ(workspace.entries.capacity(),entries_capacity);
    EXPECT_EQ(workspace.response.capacity(),response_capacity);
    EXPECT_EQ(workspace.local_column.capacity(),local_column_capacity);
    EXPECT_EQ(workspace.tile_columns.capacity(),tile_columns_capacity);
    compare_bitwise(fresh.beta,repeated.beta);
}

TEST(JointComponentNumericsTest, UnweightedLinearSolveRejectsInvalidDimensionsAndNonfiniteInputs)
{
    namespace n=p::runtime;
    const auto expect_invalid=[](const n::Sparse & design,const Vector & observations) {
        const auto result=n::SolveLinear(design,observations);
        EXPECT_FALSE(result.valid);
        EXPECT_EQ(result.reason,"invalid-input");
    };

    const n::Sparse empty_columns(6,0);
    expect_invalid(empty_columns,Vector::Ones(6));
    const Matrix odd_columns=Matrix::Zero(6,3);
    expect_invalid(n::Sparse(odd_columns.sparseView()),Vector::Ones(6));
    const Matrix too_few_rows=Matrix::Identity(4,4);
    expect_invalid(n::Sparse(too_few_rows.sparseView()),Vector::Ones(4));

    Matrix valid_design=Matrix::Zero(6,4);
    for(Eigen::Index k=0;k<4;++k) valid_design(k,k)=1;
    const n::Sparse valid_sparse(valid_design.sparseView());
    expect_invalid(valid_sparse,Vector::Ones(5));

    Matrix nonfinite_design=valid_design;
    nonfinite_design(0,0)=std::numeric_limits<double>::infinity();
    expect_invalid(n::Sparse(nonfinite_design.sparseView()),Vector::Ones(6));
    Vector nonfinite_observations=Vector::Ones(6);
    nonfinite_observations(0)=std::numeric_limits<double>::quiet_NaN();
    expect_invalid(valid_sparse,nonfinite_observations);
}

TEST(JointComponentNumericsTest, LocalLinearCertificateKeepsItsStrictKktLimit)
{
    n::Sparse design(2,2);
    design.insert(0,0)=1.0; design.insert(1,1)=1.0; design.makeCompressed();
    Vector beta=Vector::Ones(2),observations(2);
    observations<<1.0-2e-7,1.0;
    const auto certificate=n::CertifyLinear(design,observations,beta,observations.norm());
    ASSERT_TRUE(certificate.available);
    EXPECT_GT(certificate.projected_kkt,n::PracticalJointConvergenceCriteria.global_ac_kkt);
    EXPECT_FALSE(certificate.kkt_passed);
}

TEST(JointComponentNumericsTest, CompactSvdMatchesJacobiReferenceAtRankBoundaryAndAbsoluteOverride)
{
    namespace n=p::runtime;
    constexpr double relative=1e-9;
    for(double absolute:{-1.,1e-7}) for(double factor:{0.,.5,.999999,1.,1.000001,2.})
    {
        const double threshold=absolute<0 ? relative : absolute;
        Matrix x=Matrix::Identity(16,16);
        for(Eigen::Index diagonal=0;diagonal<15;++diagonal)
            x(diagonal,diagonal)=1.-.01*static_cast<double>(diagonal);
        x(15,15)=threshold*factor;
        Eigen::JacobiSVD<Matrix,Eigen::ComputeThinV> reference(x);
        const double maximum=reference.singularValues()(0);
        reference.setThreshold(absolute>=0 && maximum>0 ? absolute/maximum : relative);
        const double expected_threshold=reference.threshold()*maximum;
        const Eigen::Index expected_rank=reference.rank();
        const auto actual=n::CompactSvd(x,relative,absolute,n::CompactSvdVectors::Right);
        ASSERT_TRUE(actual.valid); EXPECT_EQ(actual.rank,expected_rank);
        EXPECT_EQ(actual.threshold,expected_threshold);
        EXPECT_LT((actual.singular_values-reference.singularValues()).norm(),1e-10);
        ASSERT_EQ(actual.right_vectors.cols(),reference.matrixV().cols());
        for(Eigen::Index column=0;column<actual.right_vectors.cols();++column)
            EXPECT_NEAR(std::abs(actual.right_vectors.col(column).dot(reference.matrixV().col(column))),1.0,1e-10);
        EXPECT_TRUE(actual.used_bdc);
        const bool near=std::abs(factor-1)*threshold<=64*std::numeric_limits<double>::epsilon()*16;
        EXPECT_EQ(actual.jacobi_retry,near);
    }
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

TEST(JointComponentNumericsTest, DenseSpectrumReferenceUsesGlobalRowsDecisionDimensionsAndStrictBoundary)
{
    namespace n=p::runtime;
    const n::RankPolicy global_rows{1000,4,2};
    constexpr Eigen::Index decisions=4;
    const double threshold=global_rows.Relative(decisions);
    Matrix boundary=Matrix::Zero(6,4);
    boundary(0,0)=1;
    boundary(1,1)=2*threshold;
    boundary(2,2)=threshold;
    boundary(3,3)=.5*threshold;

    const auto dense=p::DenseSpectrum(boundary,global_rows,decisions,false);
    ASSERT_TRUE(dense.available)<<dense.reason;
    EXPECT_EQ(dense.rows,global_rows.rows);
    EXPECT_EQ(dense.columns,decisions);
    EXPECT_EQ(dense.full_design_columns,decisions);
    EXPECT_DOUBLE_EQ(dense.threshold,threshold);
    EXPECT_EQ(dense.rank,2);
    EXPECT_EQ(dense.singular_values.size(),4);
    EXPECT_DOUBLE_EQ(dense.singular_values(0),1);
    EXPECT_DOUBLE_EQ(dense.singular_values(1),2*threshold);
    EXPECT_DOUBLE_EQ(dense.singular_values(2),threshold);
    EXPECT_DOUBLE_EQ(dense.singular_values(3),.5*threshold);

    const auto native=n::EvaluateRank(boundary,{global_rows,decisions,-1,n::RankBoundary::SvdNative});
    const auto strict=n::EvaluateRank(boundary,{global_rows,decisions,-1,n::RankBoundary::StrictGreater});
    ASSERT_TRUE(native.valid && strict.valid);
    EXPECT_EQ(native.rank,3);
    EXPECT_EQ(strict.rank,2);
    EXPECT_EQ(dense.rank,strict.rank);

    const auto normalized=p::DenseSpectrum(boundary,global_rows,decisions,true);
    ASSERT_TRUE(normalized.available)<<normalized.reason;
    EXPECT_EQ(normalized.rank,4);
    EXPECT_LE((normalized.singular_values-Vector::Ones(4)).norm(),1e-14);
    EXPECT_NEAR(normalized.column_norms(1),2*threshold,1e-27);

    Matrix rank_deficient=Matrix::Zero(6,4);
    rank_deficient.diagonal()<<1,.5,.25,0;
    const auto deficient=p::DenseSpectrum(rank_deficient,global_rows,decisions,false);
    ASSERT_TRUE(deficient.available);
    EXPECT_EQ(deficient.rank,3);
    rank_deficient(3,3)=2*threshold;
    const auto near_deficient=p::DenseSpectrum(rank_deficient,global_rows,decisions,false);
    ASSERT_TRUE(near_deficient.available);
    EXPECT_EQ(near_deficient.rank,4);

    const n::RankPolicy decision_columns{4,1000,4};
    Matrix decision_limited=Matrix::Zero(6,4);
    decision_limited.diagonal()<<1,.5,.25,100*std::numeric_limits<double>::epsilon();
    const auto decision_reference=p::DenseSpectrum(decision_limited,decision_columns,1000,false);
    ASSERT_TRUE(decision_reference.available);
    EXPECT_DOUBLE_EQ(decision_reference.threshold,1000*std::numeric_limits<double>::epsilon());
    EXPECT_EQ(decision_reference.rank,3);
}

TEST(JointComponentNumericsTest, MemoryPlannerRetainsNonRankStagesAndLocalBudgets)
{
    using rhbm_gem::core::JointProblemInput;
    std::size_t bytes{};
    EXPECT_TRUE(n::CheckedMultiplyBytes(8,16,bytes)); EXPECT_EQ(bytes,128u);
    EXPECT_FALSE(n::CheckedMultiplyBytes(std::numeric_limits<std::size_t>::max(),2,bytes));
    EXPECT_FALSE(n::CheckedAddBytes(std::numeric_limits<std::size_t>::max(),1,bytes));

    JointProblemInput input; input.observations={1,2}; input.atom_ids={"a","b"}; input.row_ids={"r0","r1"};
    input.support={{{0,0},{1,1}},{{1,2}}};
    ASSERT_EQ(n::CountSupportEntries(input),3u);
    ASSERT_TRUE(n::EstimateSnapshotBytes(input));
    ASSERT_TRUE(n::EstimateConstructionScratchBytes(input));

    n::MemoryShape large; large.component_atoms=37'406; large.component_observations=5'938'087;
    large.component_support_entries=2'000'000; large.component_blocks=3'118;
    large.affected_row_entries=3'000'000; large.component_free_columns=74'812;
    large.maximum_block_rows=8'000; large.maximum_block_support_entries=80'000;
    large.maximum_neighbor_atoms=12; large.component_width_columns=12;
    large.maximum_tile_rows=8192; large.resident_bytes=256'000'000;
    const auto construction=n::PlanMemory(large,n::MemoryStage::ProblemConstruction);
    const auto structural=n::PlanMemory(large,n::MemoryStage::StructuralConstruction);
    const auto preparation=n::PlanMemory(large,n::MemoryStage::FixedNeighborPreparation);
    const auto state=n::PlanMemory(large,n::MemoryStage::GlobalStateEvaluation);
    const auto endpoint=n::PlanMemory(large,n::MemoryStage::EndpointVerification);
    auto local=large; local.component_free_columns=24;
    const auto conditioning=n::PlanMemory(local,n::MemoryStage::LocalConditioning);
    EXPECT_EQ(construction.status,n::MemoryStatus::Pass);
    EXPECT_EQ(structural.status,n::MemoryStatus::Pass);
    EXPECT_EQ(preparation.status,n::MemoryStatus::Pass);
    EXPECT_EQ(state.status,n::MemoryStatus::Pass);
    EXPECT_EQ(endpoint.status,n::MemoryStatus::Pass);
    EXPECT_EQ(conditioning.status,n::MemoryStatus::Pass);
    EXPECT_EQ(conditioning.factor_rows,24);
    EXPECT_EQ(conditioning.factor_columns,24);
    EXPECT_EQ(conditioning.budget_bytes,n::JointMemoryPlanningBudgetBytes);

    auto bounded_affected_replay=large;
    bounded_affected_replay.affected_replay_scratch=true;
    bounded_affected_replay.maximum_block_replay_entries=200'000;
    bounded_affected_replay.maximum_block_atoms=24;
    const auto full_replay=n::PlanMemory(large,n::MemoryStage::GlobalStateEvaluation);
    const auto d1_at_full_budget=n::PlanMemory(bounded_affected_replay,n::MemoryStage::GlobalStateEvaluation,
        full_replay.estimated_peak_bytes);
    auto fallback_shape=bounded_affected_replay; fallback_shape.affected_replay_scratch=false;
    const auto fallback_at_full_budget=n::PlanMemory(fallback_shape,n::MemoryStage::GlobalStateEvaluation,
        full_replay.estimated_peak_bytes);
    EXPECT_EQ(full_replay.status,n::MemoryStatus::Pass);
    EXPECT_EQ(d1_at_full_budget.status,n::MemoryStatus::Reject);
    EXPECT_EQ(fallback_at_full_budget.status,n::MemoryStatus::Pass);
    EXPECT_GT(d1_at_full_budget.estimated_peak_bytes,fallback_at_full_budget.estimated_peak_bytes);

    auto overflow=large; overflow.component_free_columns=std::numeric_limits<std::size_t>::max();
    EXPECT_EQ(n::PlanMemory(overflow,n::MemoryStage::LocalConditioning).status,n::MemoryStatus::Overflow);
    auto tile_overflow=large; tile_overflow.maximum_tile_rows=std::numeric_limits<std::size_t>::max();
    EXPECT_EQ(n::PlanMemory(tile_overflow,n::MemoryStage::LocalConditioning).status,n::MemoryStatus::Overflow);
    auto construction_overflow=large;
    construction_overflow.component_observations=std::numeric_limits<std::size_t>::max();
    EXPECT_EQ(n::PlanMemory(construction_overflow,n::MemoryStage::ProblemConstruction).status,n::MemoryStatus::Overflow);
}

TEST(JointComponentNumericsTest, StreamingStateMatchesSparseGlobalEvaluation)
{
    using rhbm_gem::core::JointProblem;
    using rhbm_gem::core::JointProblemAccess;
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const auto & data=JointProblemAccess::Get(problem);
    const Vector eta=Vector::Constant(8,std::log(.55));
    Vector beta(16);
    for(Eigen::Index atom=0;atom<8;++atom)
    {
        beta(2*atom)=atom%5==0 ? 0.0 : 1.7+.1*static_cast<double>(atom);
        beta(2*atom+1)=atom%2 ? -.4 : .3;
    }
    const auto sparse=n::EvaluateState(data.domain,data.y,eta,beta,data.context);
    const auto streaming=n::EvaluateStateStreaming(data.domain,data.y,eta,beta,data.context);
    ASSERT_TRUE(sparse.valid)<<sparse.reason;
    ASSERT_TRUE(streaming.valid)<<streaming.reason;
    ASSERT_TRUE(sparse.certificate.available && streaming.certificate.available);
    const Vector sparse_prediction=sparse.residual+data.y;
    const Vector streaming_prediction=streaming.residual+data.y;
    const double prediction_scale=std::max(1.0,sparse_prediction.cwiseAbs().maxCoeff());
    EXPECT_LE((streaming_prediction-sparse_prediction).lpNorm<Eigen::Infinity>(),2e-13*prediction_scale);
    EXPECT_LE((streaming.residual-sparse.residual).lpNorm<Eigen::Infinity>(),2e-13*prediction_scale);
    EXPECT_NEAR(streaming.certificate.objective,sparse.certificate.objective,2e-13*std::max(1.0,sparse.certificate.objective));
    EXPECT_NEAR(streaming.certificate.projected_kkt,sparse.certificate.projected_kkt,2e-12);
    EXPECT_EQ(streaming.certificate.active_atoms,sparse.certificate.active_atoms);
    EXPECT_EQ(streaming.certificate.feasible,sparse.certificate.feasible);
    EXPECT_LE((streaming.gradient-sparse.gradient).lpNorm<Eigen::Infinity>(),
        2e-12*std::max(1.0,sparse.gradient.lpNorm<Eigen::Infinity>()));
    EXPECT_EQ(streaming.x.size(),0);
    EXPECT_EQ(streaming.derivative.size(),0);
}

TEST(JointComponentNumericsTest, IndependentStreamingReplayMatchesLegacyMatrixReplay)
{
    for(const std::string topology:{"chain","cube"})
    {
        SCOPED_TRACE(topology);
        using rhbm_gem::core::JointProblem;
        using rhbm_gem::core::JointProblemAccess;
        JointProblem problem(second_stage_test::SyntheticJointWorkload(topology,8));
        const auto & data=JointProblemAccess::Get(problem);
        const Vector eta=Vector::Constant(8,std::log(.55));
        const auto primary=n::EvaluateProfile(data.domain,data.y,eta,&data.context);
        ASSERT_TRUE(primary.valid);
        const auto streaming=n::EvaluateStateStreaming(data.domain,data.y,primary.eta,primary.beta,data.context);
        ASSERT_TRUE(streaming.valid);
        const auto legacy=n::CheckReplay(data.domain,data.y,primary,data.context);
        const auto independent=n::CheckIndependentReplay(data.domain,data.y,streaming,data.context);
        ASSERT_TRUE(legacy.passed)<<legacy.reason;
        ASSERT_TRUE(independent.available && independent.valid)<<independent.reason;
        EXPECT_TRUE(independent.primary_kkt_passed);
        EXPECT_TRUE(independent.independent_kkt_passed);
        EXPECT_TRUE(independent.prediction_passed);
        EXPECT_TRUE(independent.residual_passed);
        EXPECT_TRUE(independent.objective_passed);
        EXPECT_TRUE(independent.kkt_agrees);
        EXPECT_TRUE(independent.gradient_agrees);
        EXPECT_TRUE(independent.passed)<<independent.reason;
        EXPECT_LE(independent.prediction_difference,2e-12+2e-13*(primary.x*primary.beta).cwiseAbs().maxCoeff());
        EXPECT_LE(independent.residual_difference,2e-12+2e-13*(primary.x*primary.beta).cwiseAbs().maxCoeff());
        EXPECT_LE(independent.objective_difference,1e-12+2e-12*primary.certificate.objective);
        EXPECT_LE(independent.kkt_difference,1e-13);
        EXPECT_LE(independent.gradient_difference,1e-13+2e-9*streaming.gradient.lpNorm<Eigen::Infinity>());
        EXPECT_TRUE(std::isfinite(independent.cancellation_ratio));
    }
}

TEST(JointComponentNumericsTest, PracticalStationarityCannotHidePerCoordinateReplayDisagreement)
{
    constexpr Eigen::Index rows=4;
    std::vector<std::vector<n::Support>> support(1);
    n::Vector y(rows),eta= n::Vector::Constant(1,std::log(.55)),beta(2);
    beta<<1.5,-.1;
    for(Eigen::Index row=0;row<rows;++row)
    {
        const double square=std::array<double,4>{0.0,.36,1.44,4.0}[static_cast<std::size_t>(row)];
        support[0].push_back({row,square});
        const auto basis=n::EvaluateKernel(square,.55,2.5);
        y(row)=beta(0)*basis.gaussian+beta(1)*basis.charge;
    }
    const n::Domain domain(rows,std::move(support));
    auto context=p::MakeContext(y,1);
    context.atom_ids=n::Identities({"target-42"});
    context.row_ids=n::Identities({"row-1","row-2","row-3","row-4"});
    auto primary=n::EvaluateStateStreaming(domain,y,eta,beta,context);
    ASSERT_TRUE(primary.valid)<<primary.reason;
    auto original=n::CheckIndependentReplay(domain,y,primary,context,nullptr,nullptr,
        n::PracticalJointConvergenceCriteria);
    ASSERT_TRUE(original.passed)<<original.reason;

    primary.gradient(0)+=5e-10;
    const auto replay=n::CheckIndependentReplay(domain,y,primary,context,nullptr,nullptr,
        n::PracticalJointConvergenceCriteria);
    ASSERT_TRUE(replay.valid)<<replay.reason;
    EXPECT_TRUE(replay.primary_kkt_passed);
    EXPECT_TRUE(replay.independent_kkt_passed);
    EXPECT_TRUE(replay.prediction_passed);
    EXPECT_TRUE(replay.residual_passed);
    EXPECT_TRUE(replay.objective_passed);
    EXPECT_TRUE(replay.kkt_agrees);
    EXPECT_FALSE(replay.gradient_agrees);
    EXPECT_FALSE(replay.passed);
    EXPECT_EQ(replay.reason,"per-coordinate-replay-enclosure");
    EXPECT_GT(replay.primary_width_gradient,n::StrictJointConvergenceCriteria.global_width_gradient);
    EXPECT_LE(replay.primary_width_gradient,n::PracticalJointConvergenceCriteria.global_width_gradient);
    EXPECT_EQ(replay.maximum_gradient_difference_atom_id,"target-42");
    EXPECT_EQ(replay.gradient_enclosure_failure_atom_id,"target-42");
    EXPECT_DOUBLE_EQ(replay.gradient_enclosure_failure_primary,primary.gradient(0));
    EXPECT_TRUE(std::isfinite(replay.gradient_enclosure_failure_independent));
    EXPECT_DOUBLE_EQ(replay.gradient_enclosure_failure_difference,
        std::abs(replay.gradient_enclosure_failure_independent-primary.gradient(0)));
    EXPECT_DOUBLE_EQ(replay.gradient_difference,replay.gradient_enclosure_failure_difference);
    EXPECT_DOUBLE_EQ(replay.maximum_gradient_difference_primary,primary.gradient(0));
    EXPECT_DOUBLE_EQ(replay.maximum_gradient_difference_independent,replay.gradient_enclosure_failure_independent);
    EXPECT_DOUBLE_EQ(replay.maximum_gradient_difference_bound,replay.gradient_enclosure_failure_bound);
    EXPECT_DOUBLE_EQ(replay.maximum_gradient_difference_ratio,replay.gradient_enclosure_failure_ratio);
    EXPECT_DOUBLE_EQ(replay.primary_width_gradient,std::abs(primary.gradient(0)));
    EXPECT_DOUBLE_EQ(replay.independent_width_gradient,std::abs(replay.gradient_enclosure_failure_independent));
    EXPECT_DOUBLE_EQ(replay.gradient_enclosure_failure_bound,
        1e-13+2e-9*std::abs(primary.gradient(0)));
    EXPECT_DOUBLE_EQ(replay.gradient_enclosure_failure_ratio,
        replay.gradient_enclosure_failure_difference/replay.gradient_enclosure_failure_bound);
    EXPECT_GT(replay.gradient_enclosure_failure_ratio,1.0);
    EXPECT_DOUBLE_EQ(replay.primary_objective,primary.certificate.objective);
    EXPECT_DOUBLE_EQ(replay.primary_residual_norm,primary.residual.norm());
    EXPECT_DOUBLE_EQ(replay.independent_residual_norm,std::sqrt(2.0*replay.independent_objective));

    const auto assessment=n::VerifyJointEndpoint(domain,y,primary,replay,context,
        n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(assessment.numerical_qualified);
    EXPECT_EQ(assessment.failure,"per-coordinate-replay-enclosure");
}

TEST(JointComponentNumericsTest, IndependentStreamingReplayPreservesActiveBoundaryAndCancellation)
{
    constexpr Eigen::Index rows=33;
    std::vector<std::vector<n::Support>> support(1);
    Vector y(rows),eta=Vector::Constant(1,std::log(.55)),beta(2);
    Vector gaussian(rows),charge(rows);
    for(Eigen::Index row=0;row<rows;++row)
    {
        const double square=6.25*static_cast<double>(row)/static_cast<double>(rows-1);
        support[0].push_back({row,square});
        const auto basis=n::EvaluateKernel(square,.55,2.5);
        gaussian(row)=basis.gaussian; charge(row)=basis.charge;
    }
    n::Domain domain(rows,std::move(support));
    y=-gaussian; beta<<0,charge.dot(y)/charge.squaredNorm();
    const auto context=p::MakeContext(y,1);
    const auto active=n::EvaluateStateStreaming(domain,y,eta,beta,context);
    ASSERT_TRUE(active.valid)<<active.reason;
    ASSERT_TRUE(active.certificate.kkt_passed);
    const auto active_replay=n::CheckIndependentReplay(domain,y,active,context);
    EXPECT_TRUE(active_replay.passed)<<active_replay.reason;
    EXPECT_TRUE(active_replay.feasible);
    EXPECT_DOUBLE_EQ(beta(0),0.0);

    std::vector<std::vector<n::Support>> duplicate_support(2);
    for(Eigen::Index row=0;row<rows;++row)
    {
        const double square=6.25*static_cast<double>(row)/static_cast<double>(rows-1);
        duplicate_support[0].push_back({row,square});
        duplicate_support[1].push_back({row,square});
    }
    n::Domain duplicate_domain(rows,std::move(duplicate_support));
    const Vector zero=Vector::Zero(rows),duplicate_eta=Vector::Constant(2,std::log(.55));
    Vector cancelling_beta(4); cancelling_beta<<0,1e12,0,-1e12;
    const auto duplicate_context=p::MakeContext(zero,2);
    const auto cancelling=n::EvaluateStateStreaming(duplicate_domain,zero,duplicate_eta,cancelling_beta,duplicate_context);
    ASSERT_TRUE(cancelling.valid)<<cancelling.reason;
    const auto cancellation_replay=n::CheckIndependentReplay(duplicate_domain,zero,cancelling,duplicate_context);
    EXPECT_TRUE(cancellation_replay.passed)<<cancellation_replay.reason;
    EXPECT_GT(cancellation_replay.cancellation_ratio,1e10);
    EXPECT_DOUBLE_EQ(cancellation_replay.objective_difference,0.0);
}

TEST(JointComponentNumericsTest, BoundedReplayUsesProjectedWidthKktAndChecksFeasibility)
{
    constexpr Eigen::Index rows=17;
    std::vector<std::vector<n::Support>> support(1);
    Vector eta=Vector::Constant(1,std::log(.55));
    for(Eigen::Index row=0;row<rows;++row)
        support[0].push_back({row,6.25*static_cast<double>(row)/static_cast<double>(rows-1)});
    const n::Domain domain(rows,std::move(support));
    auto evaluate=[&](double truth_width) {
        Vector y(rows);
        for(Eigen::Index row=0;row<rows;++row)
        {
            const double square=6.25*static_cast<double>(row)/static_cast<double>(rows-1);
            const auto basis=n::EvaluateKernel(square,truth_width,2.5);
            y(row)=2.0*basis.gaussian+.5*basis.charge;
        }
        auto context=p::MakeContext(y,1);
        context.eta_upper_bound=std::log(.55);
        const auto primary=n::EvaluateProfile(domain,y,eta,&context);
        EXPECT_TRUE(primary.valid)<<primary.reason;
        const auto replay=n::CheckIndependentReplay(domain,y,primary,context,nullptr,nullptr,
            n::StrictJointConvergenceCriteria);
        return std::tuple{y,context,primary,replay};
    };

    auto [outward_y,outward_context,outward,outward_replay]=evaluate(1.2);
    ASSERT_TRUE(outward.valid);
    ASSERT_TRUE(outward_replay.available && outward_replay.valid)<<outward_replay.reason;
    EXPECT_LT(outward.gradient(0),0.0);
    EXPECT_TRUE(outward_replay.width_bound_applied);
    EXPECT_TRUE(outward_replay.width_bound_feasible);
    EXPECT_TRUE(outward_replay.primary_width_kkt_passed);
    EXPECT_TRUE(outward_replay.independent_width_kkt_passed)<<outward_replay.reason;
    EXPECT_GT(outward_replay.primary_width_gradient,1e-12);
    EXPECT_DOUBLE_EQ(outward_replay.primary_width_projected_kkt,0.0);
    EXPECT_DOUBLE_EQ(outward_replay.independent_width_projected_kkt,0.0);
    EXPECT_DOUBLE_EQ(*n::ProfileGradientInfinityNorm(outward,outward_context),0.0);
    EXPECT_EQ(outward_replay.active_width_bounds,1u);
    EXPECT_GT(outward_replay.boundary_multiplier_max,0.0);
    EXPECT_LE(outward_replay.boundary_gradient_max,0.0);
    EXPECT_DOUBLE_EQ(outward_replay.complementarity_residual,0.0);
    const auto bounded_assessment=n::VerifyJointEndpoint(domain,outward_y,outward,outward_replay,outward_context);
    EXPECT_TRUE(bounded_assessment.numerical_qualified)<<bounded_assessment.failure;
    const auto bounded_evidence=n::AssessmentEvidence(bounded_assessment,
        rhbm_gem::JointEvidenceScope::AssembledGlobal);
    const auto width_stationarity=std::find_if(bounded_evidence.begin(),bounded_evidence.end(),[](const auto & check) {
        return check.name=="width-stationarity";
    });
    ASSERT_NE(width_stationarity,bounded_evidence.end());
    EXPECT_EQ(width_stationarity->status,rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(width_stationarity->value.has_value());
    EXPECT_DOUBLE_EQ(*width_stationarity->value,0.0);
    const auto bound_feasibility=std::find_if(bounded_evidence.begin(),bounded_evidence.end(),[](const auto & check) {
        return check.name=="width-bound-feasibility";
    });
    ASSERT_NE(bound_feasibility,bounded_evidence.end());
    EXPECT_EQ(bound_feasibility->status,rhbm_gem::JointCheckStatus::Passed);

    auto [inward_y,inward_context,inward,inward_replay]=evaluate(.3);
    ASSERT_TRUE(inward.valid);
    ASSERT_TRUE(inward_replay.available && inward_replay.valid)<<inward_replay.reason;
    EXPECT_GT(inward.gradient(0),0.0);
    EXPECT_FALSE(inward_replay.primary_width_kkt_passed);
    EXPECT_FALSE(inward_replay.independent_width_kkt_passed);
    EXPECT_FALSE(inward_replay.passed);
    const auto rejected=n::VerifyJointEndpoint(domain,inward_y,inward,inward_replay,inward_context);
    EXPECT_FALSE(rejected.numerical_qualified);

    const Vector infeasible_eta=Vector::Constant(1,std::log(.55)+1e-8);
    const auto infeasible=n::EvaluateProfile(domain,outward_y,infeasible_eta,&outward_context);
    ASSERT_TRUE(infeasible.valid);
    const auto infeasible_replay=n::CheckIndependentReplay(domain,outward_y,infeasible,outward_context);
    ASSERT_TRUE(infeasible_replay.available && infeasible_replay.valid);
    EXPECT_FALSE(infeasible_replay.width_bound_feasible);
    EXPECT_FALSE(infeasible_replay.passed);
    EXPECT_NEAR(infeasible_replay.maximum_width_bound_violation,1e-8,1e-12);
    const auto infeasible_assessment=n::VerifyJointEndpoint(domain,outward_y,infeasible,
        infeasible_replay,outward_context);
    const auto infeasible_evidence=n::AssessmentEvidence(infeasible_assessment,
        rhbm_gem::JointEvidenceScope::AssembledGlobal);
    const auto infeasible_bound=std::find_if(infeasible_evidence.begin(),infeasible_evidence.end(),[](const auto & check) {
        return check.name=="width-bound-feasibility";
    });
    ASSERT_NE(infeasible_bound,infeasible_evidence.end());
    EXPECT_EQ(infeasible_bound->status,rhbm_gem::JointCheckStatus::Failed);
    const auto infeasible_trial_trust=n::CheckReplay(domain,outward_y,infeasible,outward_context);
    EXPECT_FALSE(infeasible_trial_trust.passed);
    EXPECT_EQ(infeasible_trial_trust.reason,"width-bound-infeasible");
}

TEST(JointComponentNumericsTest, WidthBoundKktUsesUpperFaceSignAndInteriorGradient)
{
    Vector eta(2),gradient(2);
    eta<<0.0,1.0;
    gradient<<0.0,-2.0;
    auto kkt=n::CalculateWidthBoundKkt(eta,gradient,1.0);
    EXPECT_TRUE(kkt.feasible);
    EXPECT_DOUBLE_EQ(kkt.projected,0.0);
    EXPECT_EQ(kkt.active,1u);
    EXPECT_DOUBLE_EQ(kkt.boundary_gradient,-2.0);
    EXPECT_DOUBLE_EQ(kkt.multiplier,2.0);
    EXPECT_DOUBLE_EQ(kkt.complementarity,0.0);

    gradient(1)=.25;
    kkt=n::CalculateWidthBoundKkt(eta,gradient,1.0);
    EXPECT_TRUE(kkt.feasible);
    EXPECT_DOUBLE_EQ(kkt.projected,.25);

    eta(0)=std::nextafter(1.0,0.0);
    eta(1)=1.0;
    gradient<<-.5,0.0;
    kkt=n::CalculateWidthBoundKkt(eta,gradient,1.0);
    EXPECT_TRUE(kkt.feasible);
    EXPECT_NEAR(kkt.projected,0.0,2e-16);

    eta(1)=1.0+1e-8;
    kkt=n::CalculateWidthBoundKkt(eta,gradient,1.0);
    EXPECT_FALSE(kkt.feasible);
    EXPECT_NEAR(kkt.violation,1e-8,1e-12);

    eta<<.2,.4;
    gradient<<-.3,.5;
    kkt=n::CalculateWidthBoundKkt(eta,gradient,1.0);
    EXPECT_TRUE(kkt.feasible);
    EXPECT_EQ(kkt.active,0u);
    EXPECT_DOUBLE_EQ(kkt.projected,gradient.lpNorm<Eigen::Infinity>());
}
