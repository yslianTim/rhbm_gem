#pragma once

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <tuple>
#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

class AtomObject;

class GroupPotentialEntryBase
{
public:
    virtual ~GroupPotentialEntryBase() = default;
    virtual void AddAtomObjectPtr(const void * group_key, AtomObject * atom_object) = 0;
    virtual void AddGausEstimateMean(const void * group_key, Eigen::VectorXd vec) = 0;
    virtual void AddGausEstimateMDPDE(const void * group_key, Eigen::VectorXd vec) = 0;
    virtual void AddGausEstimatePrior(const void * group_key, Eigen::VectorXd vec) = 0;
    virtual void AddGausVariancePrior(const void * group_key, Eigen::VectorXd vec) = 0;
    virtual const std::vector<AtomObject *> & GetAtomObjectPtrList(const void * group_key) const = 0;
    virtual std::tuple<double, double> GetGausEstimateMean(const void * group_key) const = 0;
    virtual std::tuple<double, double> GetGausEstimateMDPDE(const void * group_key) const = 0;
    virtual std::tuple<double, double> GetGausEstimatePrior(const void * group_key) const = 0;
    virtual std::tuple<double, double> GetGausVariancePrior(const void * group_key) const = 0;
};