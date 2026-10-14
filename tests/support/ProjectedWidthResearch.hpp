#pragma once
#include "core/detail/joint_component/Numerics.hpp"

namespace second_stage_test::projected_width {
namespace runtime=rhbm_gem::core::joint_component;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;

struct Fixture
{
    std::string name;
    Matrix free_design,width_derivative;
    runtime::RankPolicy policy;
    Eigen::Index full_design_columns{},target_columns{};
    bool active_amplitude_boundary{};
};

struct SpectrumOracle
{
    Vector singular_values;
    Eigen::Index rank{};
    double threshold{},minimum{};
    bool threshold_sensitive{},full_rank{};
};

struct DenseOracle
{
    SpectrumOracle design,width;
    Vector design_column_norms;
    Matrix projected_width;
    Eigen::Index rows{},free_columns{},width_columns{},structural_nonzeros{},supported_rows{};
    double projection_residual{};
    bool valid{};
    std::string reason;
};

std::vector<Fixture> Fixtures();
DenseOracle EvaluateDenseOracle(const Fixture &);
}
