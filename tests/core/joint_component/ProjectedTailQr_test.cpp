#include <gtest/gtest.h>
#include "core/detail/joint_component/CompactSvd.hpp"
#include "core/detail/joint_component/SparseFactor.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {
namespace n=rhbm_gem::core::joint_component;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;
using Sparse=n::Sparse;

struct TailFixture
{
    std::string name;
    Matrix design,raw,residual,correction;
    std::vector<Eigen::Index> free_columns,width_order;
    double scale{1.};
};

Sparse SparseOf(const Matrix & value)
{return value.sparseView(0.,0.);}

Matrix ReorderColumns(const Matrix & value,const std::vector<Eigen::Index> & order)
{
    Matrix result(value.rows(),static_cast<Eigen::Index>(order.size()));
    for(Eigen::Index k=0;k<result.cols();++k) result.col(k)=value.col(order[static_cast<std::size_t>(k)]);
    return result;
}

TailFixture MakeFixture(const std::string & name,bool cube=false,bool active_face=false,
    bool near_collinear=false,double cancellation=1.,bool row_permutation=false,
    bool column_scaling=false,bool near_rank_boundary=false,bool large_dynamic_range=false,
    bool permute_widths=false)
{
    constexpr Eigen::Index rows=56,m=4;
    const Eigen::Index p=active_face ? 6 : 8;
    Matrix z=Matrix::Zero(rows,p),width=Matrix::Zero(rows,m);
    for(Eigen::Index column=0;column<p;++column)
    {
        z(column,column)=1.+.03*static_cast<double>(column);
        if(cube)
        {
            z(p+(column+1)%p,column)=.17;
            z(2*p+column,column)=.23;
            z(3*p+(3*column+1)%p,column)=.08;
        }
        else
        {
            if(column+1<p) z(column+1,column)=.19;
            z(p+column,column)=.14;
            z(2*p+column,column)=.06;
        }
    }
    for(Eigen::Index column=0;column<m;++column)
    {
        width(32+column,column)=1.+.2*static_cast<double>(column);
        width(40+(column+1)%m,column)=.25;
        width(48+(2*column+1)%m,column)=.11;
    }
    if(near_collinear)
    {
        const Matrix original=width;
        width.col(1)=original.col(0)+1e-6*original.col(1);
    }
    if(near_rank_boundary)
    {
        const Matrix original=width;
        const double cutoff=std::numeric_limits<double>::epsilon()*static_cast<double>(rows);
        width.col(1)=original.col(0)+.25*cutoff*original.col(1);
    }
    width*=cancellation;
    Matrix coefficients(p,m);
    for(Eigen::Index row=0;row<p;++row) for(Eigen::Index column=0;column<m;++column)
        coefficients(row,column)=std::sin(.17*static_cast<double>(row+1)*static_cast<double>(column+2));
    Matrix raw=z*coefficients+width;
    if(column_scaling || large_dynamic_range)
    {
        const std::array<double,4> moderate{{.01,.1,10.,100.}};
        const std::array<double,4> wide{{1e-6,1e-2,1e2,1e6}};
        for(Eigen::Index column=0;column<m;++column)
            raw.col(column)*=(large_dynamic_range ? wide[static_cast<std::size_t>(column)] :
                moderate[static_cast<std::size_t>(column)]);
    }
    Matrix residual(rows,1),correction(p,m);
    for(Eigen::Index row=0;row<rows;++row)
        residual(row)=std::cos(.13*static_cast<double>(row))+.01*static_cast<double>(row);
    for(Eigen::Index row=0;row<p;++row) for(Eigen::Index column=0;column<m;++column)
        correction(row,column)=.03*std::cos(.11*static_cast<double>(row+2*column));
    if(row_permutation)
    {
        Matrix permuted_z(rows,p),permuted_raw(rows,m),permuted_residual(rows,1);
        for(Eigen::Index row=0;row<rows;++row)
        {
            const Eigen::Index source=(row*17)%rows;
            permuted_z.row(row)=z.row(source); permuted_raw.row(row)=raw.row(source);
            permuted_residual.row(row)=residual.row(source);
        }
        z=std::move(permuted_z); raw=std::move(permuted_raw); residual=std::move(permuted_residual);
    }
    std::vector<Eigen::Index> free_columns;
    if(active_face) free_columns={1,2,3,5,6,7};
    else for(Eigen::Index column=0;column<p;++column) free_columns.push_back(column);
    std::vector<Eigen::Index> width_order{0,1,2,3};
    if(permute_widths) width_order={2,0,3,1};
    return {name,std::move(z),std::move(raw),std::move(residual),std::move(correction),
        std::move(free_columns),std::move(width_order),1.7};
}

Matrix RestoreWidthColumns(const Matrix & permuted,const std::vector<Eigen::Index> & order)
{
    Matrix result=Matrix::Zero(permuted.rows(),permuted.cols());
    for(Eigen::Index column=0;column<permuted.cols();++column)
        result.col(order[static_cast<std::size_t>(column)])=permuted.col(column);
    return result;
}

void ExpectRelativeMatrixNear(const Matrix & actual,const Matrix & expected,double tolerance)
{
    ASSERT_EQ(actual.rows(),expected.rows()); ASSERT_EQ(actual.cols(),expected.cols());
    EXPECT_LE((actual-expected).norm(),tolerance*(1.+expected.norm()));
}

void ValidateTailFixture(const TailFixture & fixture)
{
    const Eigen::Index nrows=fixture.design.rows(),p=fixture.design.cols(),m=fixture.raw.cols();
    const auto z=SparseOf(fixture.design),raw=SparseOf(fixture.raw);
    const auto design_factor=n::FreeDesignFactor::Fixed(z,fixture.free_columns);
    ASSERT_EQ(design_factor->Rank(),p);
    const auto transformed=design_factor->OrthogonalTransposeTailSparseForTesting(raw,true);
    ASSERT_EQ(transformed.tail.rows(),nrows-p); ASSERT_EQ(transformed.tail.cols(),m);
    EXPECT_EQ(transformed.tail.nonZeros(),static_cast<Eigen::Index>(transformed.tail_nonzeros));

    const Matrix q_transformed=design_factor->OrthogonalTransposeForTesting(fixture.raw);
    const Matrix qz=design_factor->OrthogonalTransposeForTesting(Matrix::Identity(nrows,nrows)).transpose();
    const Matrix tail=q_transformed.bottomRows(nrows-p);
    ExpectRelativeMatrixNear(Matrix(transformed.tail),tail,3e-12);
    const Matrix q2=qz.rightCols(nrows-p);
    const Matrix projected=design_factor->ProjectComplement(fixture.raw)/fixture.scale;
    ExpectRelativeMatrixNear(projected,q2*tail/fixture.scale,3e-11);

    const auto tail_input=ReorderColumns(Matrix(transformed.tail),fixture.width_order);
    std::vector<Eigen::Index> tail_columns;
    for(const auto original:fixture.width_order) tail_columns.push_back(original);
    const auto tail_factor=n::FreeDesignFactor::Fixed(SparseOf(tail_input),tail_columns);
    ASSERT_EQ(tail_factor->Rank(),m);
    if(fixture.name=="width-column-permutation")
    {
        const auto view=tail_factor->RankView();
        ASSERT_TRUE(view);
        const bool permuted=std::any_of(view->permutation.begin(),view->permutation.end(),
            [index=std::size_t{}](const int64_t value) mutable {return value!=static_cast<int64_t>(index++);});
        EXPECT_TRUE(permuted);
    }
    const Matrix factor_permuted=tail_factor->Compact()/fixture.scale;
    const Matrix factor=RestoreWidthColumns(factor_permuted,fixture.width_order);
    const Matrix projected_gram=projected.transpose()*projected;
    const Matrix factor_gram=factor.transpose()*factor;
    ExpectRelativeMatrixNear(projected_gram,factor_gram,2e-9);

    const n::RankRequest request{{nrows,2*m,m},m,-1.,n::RankBoundary::StrictGreater};
    const auto projected_svd=n::EvaluateRank(projected,request,nullptr,n::CompactSvdVectors::Right);
    const auto factor_svd=n::EvaluateRank(factor,request,nullptr,n::CompactSvdVectors::Right);
    ASSERT_TRUE(projected_svd.valid); ASSERT_TRUE(factor_svd.valid);
    EXPECT_EQ(projected_svd.rank,factor_svd.rank);
    ExpectRelativeMatrixNear(projected_svd.singular_values,factor_svd.singular_values,2e-9);
    const auto projected_norms=projected.colwise().norm().transpose();
    const auto factor_norms=factor.colwise().norm().transpose();
    ExpectRelativeMatrixNear(projected_norms,factor_norms,2e-9);
    Matrix normalized_projected=projected,normalized_factor=factor;
    for(Eigen::Index column=0;column<m;++column)
    {
        ASSERT_GT(projected_norms(column),0.); ASSERT_GT(factor_norms(column),0.);
        normalized_projected.col(column)/=projected_norms(column);
        normalized_factor.col(column)/=factor_norms(column);
    }
    const auto normalized_projected_svd=n::EvaluateRank(normalized_projected,request,nullptr,
        n::CompactSvdVectors::Right);
    const auto normalized_factor_svd=n::EvaluateRank(normalized_factor,request,nullptr,
        n::CompactSvdVectors::Right);
    ASSERT_TRUE(normalized_projected_svd.valid); ASSERT_TRUE(normalized_factor_svd.valid);
    EXPECT_EQ(normalized_projected_svd.rank,normalized_factor_svd.rank);
    ExpectRelativeMatrixNear(normalized_projected_svd.singular_values,
        normalized_factor_svd.singular_values,2e-9);
    const Matrix projected_weak=projected_svd.right_vectors.rightCols(2)*
        projected_svd.right_vectors.rightCols(2).transpose();
    const Matrix factor_weak=factor_svd.right_vectors.rightCols(2)*
        factor_svd.right_vectors.rightCols(2).transpose();
    ExpectRelativeMatrixNear(projected_weak,factor_weak,2e-7);

    const Matrix transformed_residual=design_factor->OrthogonalTransposeForTesting(
        Matrix(fixture.residual/fixture.scale));
    const Matrix tail_rhs=transformed_residual.bottomRows(nrows-p);
    const Matrix q_permuted=tail_factor->OrthogonalTransposeForTesting(tail_rhs).topRows(m);
    const Vector q=q_permuted.col(0);
    const Vector projected_response=projected.transpose()*(fixture.residual/fixture.scale);
    const Vector tail_response=factor.transpose()*q;
    EXPECT_LE((projected_response-tail_response).norm(),2e-9*(1.+projected_response.norm()));

    const Matrix compact_design=design_factor->Compact();
    const Matrix jacobian=projected-fixture.design*(fixture.correction/fixture.scale);
    Matrix compact_jacobian(m+p,m);
    compact_jacobian.topRows(m)=factor;
    compact_jacobian.bottomRows(p)=-(compact_design*(fixture.correction/fixture.scale));
    Vector compact_response(m+p);
    compact_response.head(m)=q;
    compact_response.tail(p)=(compact_design*design_factor->LeastSquares(
        Matrix(fixture.residual/fixture.scale))).col(0);
    ExpectRelativeMatrixNear(jacobian.transpose()*jacobian,
        compact_jacobian.transpose()*compact_jacobian,3e-9);
    const auto jacobian_svd=n::EvaluateRank(jacobian,request,nullptr,n::CompactSvdVectors::Right);
    const auto compact_jacobian_svd=n::EvaluateRank(compact_jacobian,request,nullptr,
        n::CompactSvdVectors::Right);
    ASSERT_TRUE(jacobian_svd.valid); ASSERT_TRUE(compact_jacobian_svd.valid);
    EXPECT_EQ(jacobian_svd.rank,compact_jacobian_svd.rank);
    ExpectRelativeMatrixNear(jacobian_svd.singular_values,
        compact_jacobian_svd.singular_values,3e-9);
    if(projected_svd.rank==m && jacobian_svd.rank==m && fixture.scale==1.7 &&
        fixture.name!="near-complete-cancellation" && fixture.name!="large-dynamic-range")
    {
        const Vector direct_response=fixture.residual/fixture.scale;
        const auto direct_correction=n::EvaluateRank(jacobian,request,
            &direct_response,n::CompactSvdVectors::None);
        const auto tail_correction=n::EvaluateRank(compact_jacobian,request,&compact_response,
            n::CompactSvdVectors::None);
        ASSERT_TRUE(direct_correction.valid); ASSERT_TRUE(tail_correction.valid);
        EXPECT_LE((direct_correction.solution-tail_correction.solution).norm(),
            2e-7*(1.+direct_correction.solution.norm()));
    }
}

TEST(JointProjectedTailTest, TailQrPreservesProjectedAlgebraAndResponse)
{
    const std::vector<TailFixture> fixtures{
        MakeFixture("small-chain"),
        MakeFixture("small-cube",true),
        MakeFixture("near-collinear",true,false,true),
        MakeFixture("near-rank-boundary",false,false,false,1.,false,false,true),
        MakeFixture("active-face",true,true),
        MakeFixture("near-complete-cancellation",false,false,false,1e-8),
        MakeFixture("row-permutation",true,false,false,1.,true),
        MakeFixture("column-scaling",false,false,false,1.,false,true),
        MakeFixture("large-dynamic-range",true,false,false,1.,false,false,false,true),
        MakeFixture("width-column-permutation",true,false,false,1.,false,false,false,false,true),
    };
    for(const auto & fixture:fixtures)
    {
        SCOPED_TRACE(fixture.name);
        ValidateTailFixture(fixture);
    }
}
}
