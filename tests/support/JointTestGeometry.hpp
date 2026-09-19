#pragma once

#include <Eigen/Dense>
#include <boost/json.hpp>
#include <array>
#include <string>
#include <vector>

namespace rhbm_gem { class MapObject; }
namespace second_stage_test::matched {

using Position = std::array<double, 3>;
struct Atom { Position position; double amplitude, width, charge; };
struct Slot { std::size_t index; Position position; double coefficient; };
struct Stencil { std::array<Slot, 64> slots; bool boundary{}; };
struct DistanceWeight { double square, coefficient; };
using Design = std::vector<std::vector<DistanceWeight>>;
struct Basis { double gaussian{}, charge{}, gaussian_log_width{}, charge_log_width{}; };

double SquareDistance(const Position &, const Position &);
Stencil MakeStencil(const rhbm_gem::MapObject & generation,
    const rhbm_gem::MapObject & sampling, const Position &);
Basis EvaluateBasis(double square, double width, double cutoff);
Design MakeDesign(const std::vector<Stencil> &, const Position &);
Eigen::MatrixXd EvaluateDesign(const Design &, double width, double cutoff);
double Predict(const Stencil &, const std::vector<Atom> &, double cutoff,
    bool quantized = false, double * absolute_sum = nullptr, double * quantization_bound = nullptr);

namespace unique_grid {
struct Voxel {std::size_t index{}; Position position{};};
struct Grid {std::vector<Voxel> voxels;};
double Direct(const Position &,const std::vector<Atom> &,double cutoff);
}
} // namespace second_stage_test::matched
