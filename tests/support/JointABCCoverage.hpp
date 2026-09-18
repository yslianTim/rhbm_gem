#pragma once
#include "support/JointABCProfile.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace second_stage_test::matched::coverage {
struct Initialization
{
    Eigen::VectorXd b;
    boost::json::object evidence;
    bool valid{};
};
double ElementWidth(int element);
std::vector<Atom> SyntheticAtoms(const std::string & name);
boost::json::object Identity(const rhbm_gem::AtomObject &);
Initialization Initialize(rhbm_gem::ModelObject &, rhbm_gem::MapObject &,
    const boost::json::array & identities);
Eigen::VectorXd InitialWidths(const Eigen::VectorXd &, const boost::json::array &, const std::string &);
void Run(const std::string & model, const std::string & map,
    const std::string & manifest, const std::string & output);
} // namespace second_stage_test::matched::coverage
