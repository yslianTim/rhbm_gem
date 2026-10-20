#include "core/detail/joint_component/Numerics.hpp"

#include <Eigen/Dense>
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace {
using rhbm_gem::GaussianModel3D;
namespace joint = rhbm_gem::core::joint_component;

struct Donor
{
    std::size_t atom_id{};
    int group{};
    std::array<double,3> position{};
    GaussianModel3D point;
    const std::array<double,3> & GetPosition() const { return position; }
};

struct Halo
{
    std::size_t atom_id{};
    int group{};
    std::array<double,3> position{};
    GaussianModel3D truth;
};

struct Scenario
{
    std::string name;
    std::vector<Donor> donors;
    std::vector<Halo> halos;
    std::vector<std::array<double,3>> raw_map_positions;
    std::vector<double> raw_map_values;
    bool clip_to_unit_box{};
};

struct Policy
{
    std::string name;
    bool group_key{};
    std::size_t k{};
    double radius{};
    bool global_fallback_on_shortage{};
};

struct Estimate
{
    std::optional<GaussianModel3D> point;
    std::size_t donor_count{};
    bool fallback{};
};

struct Summary
{
    std::size_t halos{}, available{}, fallback_count{}, donor_total{}, invalid_estimates{}, unavailable_nodes{};
    double signal_square_sum{}, parameter_a_square_sum{}, parameter_b_square_sum{},
        parameter_c_square_sum{};
    double cpu_ms{};
    std::size_t signal_count{}, observation_design_rank{};
    std::optional<double> per_atom_signal_rmse_median,per_atom_signal_rmse_max;
};

struct ObservationFit
{
    std::vector<Estimate> estimates;
    std::size_t design_rank{};
};

double SquaredDistance(const std::array<double,3> & a, const std::array<double,3> & b)
{
    double out{};
    for (std::size_t axis=0;axis<3;++axis)
    {
        const auto delta=a[axis]-b[axis];
        out+=delta*delta;
    }
    return out;
}

double Contribution(const GaussianModel3D & point, const std::array<double,3> & center,
    const std::array<double,3> & position)
{
    const auto basis=joint::EvaluateKernel(SquaredDistance(position,center),point.GetWidth(),2.5);
    return point.GetAmplitude()*basis.gaussian+point.GetOffset()*basis.charge;
}

void BuildRawMap(Scenario & scenario)
{
    std::set<std::array<double,3>> nodes;
    for (const auto & halo:scenario.halos)
        for (int ix=-5;ix<=5;++ix) for (int iy=-5;iy<=5;++iy) for (int iz=-5;iz<=5;++iz)
        {
            const std::array<double,3> position{
                halo.position[0]+0.5*ix,halo.position[1]+0.5*iy,halo.position[2]+0.5*iz};
            if (SquaredDistance(position,halo.position)>6.25) continue;
            if (scenario.clip_to_unit_box &&
                std::any_of(position.begin(),position.end(),[](double v){return v<0.0 || v>10.0;})) continue;
            nodes.insert(position);
        }
    scenario.raw_map_positions.assign(nodes.begin(),nodes.end());
    scenario.raw_map_values.reserve(scenario.raw_map_positions.size());
    for (const auto & position:scenario.raw_map_positions)
    {
        double value{};
        for (const auto & donor:scenario.donors)
            value+=Contribution(donor.point,donor.position,position);
        for (const auto & halo:scenario.halos)
            value+=Contribution(halo.truth,halo.position,position);
        scenario.raw_map_values.push_back(value);
    }
}

GaussianModel3D Median(const std::vector<const Donor *> & donors)
{
    if (donors.empty()) throw std::invalid_argument("empty donor set");
    const auto component_median=[&](auto getter) {
        std::vector<double> values;
        values.reserve(donors.size());
        for (const auto * donor:donors) values.push_back(getter(donor->point));
        std::sort(values.begin(),values.end());
        const auto middle=values.size()/2;
        return values.size()%2 ? values[middle] : (values[middle-1]+values[middle])/2.0;
    };
    return {component_median([](const auto & p){return p.GetAmplitude();}),
        component_median([](const auto & p){return p.GetWidth();}),
        component_median([](const auto & p){return p.GetOffset();})};
}

std::vector<Policy> Policies()
{
    return {
        {"group-key-global-fallback",true,0,0,false},
        {"global-target-median",false,0,0,false},
        {"knn-1-unbounded",false,1,0,false},
        {"knn-2-unbounded",false,2,0,false},
        {"knn-4-unbounded",false,4,0,false},
        {"knn-8-unbounded",false,8,0,false},
        {"knn-2-r2.5-use-available",false,2,2.5,false},
        {"knn-4-r2.5-use-available",false,4,2.5,false},
        {"knn-2-r2.5-global-shortage",false,2,2.5,true},
        {"knn-4-r5-use-available",false,4,5.0,false},
        {"knn-4-r5-global-shortage",false,4,5.0,true}
    };
}

Scenario MakeScenario(std::size_t scenario)
{
    static const std::array<std::string,10> names{{
        "nearby-similar", "halo-target-mismatch", "chemical-group-misleading",
        "same-truth-across-groups", "dense-local-cluster", "sparse-donors",
        "overlapping-halos", "map-boundary", "mixed-charge-sign", "halo-width-shift"}};
    Scenario out;
    out.name=names.at(scenario);
    const std::size_t count=scenario==4 ? 24 : scenario==5 ? 4 : 12;
    for (std::size_t i=0;i<count;++i)
    {
        const bool left=i<count/2;
        const auto cluster_index=left ? i : i-count/2;
        double x=left ? 0.25+0.24*static_cast<double>(cluster_index%6) :
            8.1+0.24*static_cast<double>(cluster_index%6);
        double y=0.25*static_cast<double>(static_cast<int>(cluster_index%3)-1);
        double z=0.2*static_cast<double>(static_cast<int>((cluster_index/3)%3)-1);
        if (scenario==4 && i<16)
        {
            x=0.35+0.10*static_cast<double>(i%4);
            y=0.10*static_cast<double>(static_cast<int>(i%3)-1);
            z=0.10*static_cast<double>(static_cast<int>((i/3)%3)-1);
        }
        if (scenario==5)
        {
            static constexpr std::array<double,4> sparse_x{{-8.0,-4.0,8.0,12.0}};
            x=sparse_x[i]; y=0.0; z=0.0;
        }
        if (scenario==7)
        {
            x=left ? 0.35+0.18*static_cast<double>(cluster_index%6) :
                9.15+0.12*static_cast<double>(cluster_index%6);
            y=left ? 0.3 : 9.4; z=left ? 0.2 : 9.6;
            out.clip_to_unit_box=true;
        }
        int group=left ? 1 : 2;
        if (scenario==2) group=left ? 2 : 1;
        double amplitude=left ? 2.0+0.12*static_cast<double>(cluster_index) :
            6.0+0.28*static_cast<double>(cluster_index);
        double width=left ? 0.55+0.02*static_cast<double>(cluster_index) :
            0.82+0.025*static_cast<double>(cluster_index);
        double charge=left ? -0.8+0.12*static_cast<double>(cluster_index) :
            1.2+0.18*static_cast<double>(cluster_index);
        if (scenario==3) {amplitude=3.0; width=0.7; charge=0.25;}
        if (scenario==4 && i<16) {amplitude=2.2; width=0.6; charge=-0.4;}
        if (scenario==8) charge=left ? -1.0-0.1*static_cast<double>(cluster_index) :
            1.0+0.1*static_cast<double>(cluster_index);
        if (scenario==9) width=0.5;
        out.donors.push_back({i+1,group,{x,y,z},{amplitude,width,charge}});
    }

    std::array<double,3> left_position{0.55,0.0,0.0};
    std::array<double,3> right_position{8.4,0.0,0.0};
    if (scenario==4) left_position={0.45,0.0,0.0};
    if (scenario==5) {left_position={0.0,0.0,0.0}; right_position={10.0,0.0,0.0};}
    if (scenario==6) {left_position={2.0,0.0,0.0}; right_position={2.0,0.0,0.0};}
    if (scenario==7) {left_position={0.1,0.1,0.1}; right_position={9.8,9.8,9.8};}
    const auto local_median=[&](int group) {
        std::vector<const Donor *> donors;
        for (const auto & donor:out.donors)
            if ((scenario==2 ? (group==1 ? donor.position[0]>4.0 : donor.position[0]<4.0) :
                    donor.group==group)) donors.push_back(&donor);
        return Median(donors);
    };
    GaussianModel3D left_truth=local_median(scenario==2 ? 2 : 1);
    GaussianModel3D right_truth=local_median(scenario==2 ? 1 : 2);
    if (scenario==1) {left_truth={15.0,0.35,-5.0}; right_truth={0.6,1.4,6.0};}
    if (scenario==2) {left_truth={2.3,0.62,-0.3}; right_truth={6.4,0.88,1.5};}
    if (scenario==3) {left_truth={3.0,0.7,0.25}; right_truth=left_truth;}
    if (scenario==4) {left_truth={2.2,0.6,-0.4}; right_truth={6.7,0.89,1.65};}
    if (scenario==5) {left_truth={4.0,0.7,0.0}; right_truth={5.0,0.8,0.5};}
    if (scenario==6) {left_truth={3.0,0.6,-0.8}; right_truth={7.5,1.05,2.0};}
    if (scenario==7) {left_truth={2.25,0.61,-0.25}; right_truth={6.9,0.9,1.7};}
    if (scenario==8) {left_truth={2.4,0.62,-1.5}; right_truth={6.6,0.9,1.6};}
    if (scenario==9) {left_truth={2.3,1.4,-0.3}; right_truth={6.5,1.5,1.5};}
    const int left_group=1;
    const int right_group=2;
    out.halos.push_back({100,left_group,left_position,left_truth});
    out.halos.push_back({101,right_group,right_position,right_truth});
    BuildRawMap(out);
    return out;
}

Estimate EstimateHalo(const Policy & policy, const Halo & halo, std::vector<Donor> & donors,
    const KDNode<Donor> * tree)
{
    std::vector<const Donor *> global;
    global.reserve(donors.size());
    for (const auto & donor:donors) global.push_back(&donor);
    const auto global_point=Median(global);
    if (policy.name=="global-target-median") return {global_point,global.size(),false};
    if (policy.group_key)
    {
        std::vector<const Donor *> group;
        for (const auto & donor:donors) if (donor.group==halo.group) group.push_back(&donor);
        if (!group.empty()) return {Median(group),group.size(),false};
        return {global_point,global.size(),true};
    }

    Donor query{halo.atom_id,halo.group,halo.position,{0.0,1.0,0.0}};
    double radius=policy.radius;
    if (radius==0.0)
    {
        const auto nearest=KDTreeAlgorithm<Donor>::KNearestNeighbors(tree,&query,policy.k);
        if (nearest.empty()) return {};
        for (const auto * donor:nearest)
            radius=std::max(radius,std::sqrt(SquaredDistance(halo.position,donor->position)));
        radius+=1e-12;
    }
    auto candidates=KDTreeAlgorithm<Donor>::RangeSearch(tree,&query,radius);
    std::sort(candidates.begin(),candidates.end(),[&](const Donor * a,const Donor * b) {
        const auto da=SquaredDistance(halo.position,a->position);
        const auto db=SquaredDistance(halo.position,b->position);
        return da==db ? a->atom_id<b->atom_id : da<db;
    });
    if (candidates.size()>policy.k) candidates.resize(policy.k);
    if (policy.global_fallback_on_shortage && candidates.size()<policy.k)
        return {global_point,global.size(),true};
    if (candidates.empty()) return {};
    const std::vector<const Donor *> selected(candidates.begin(),candidates.end());
    return {Median(selected),selected.size(),false};
}

double FixedTargetSignal(const Scenario & scenario, const std::array<double,3> & position)
{
    std::vector<const Donor *> donors;
    donors.reserve(scenario.donors.size());
    for (const auto & donor:scenario.donors) donors.push_back(&donor);
    std::sort(donors.begin(),donors.end(),[](const Donor * a,const Donor * b){return a->atom_id<b->atom_id;});
    double signal{};
    for (const auto * donor:donors) signal+=Contribution(donor->point,donor->position,position);
    return signal;
}

Summary ScoreEstimates(const Scenario & scenario, const std::vector<Estimate> & estimates,
    double cpu_ms, std::size_t design_rank=0)
{
    Summary out;
    out.halos=scenario.halos.size();
    out.cpu_ms=cpu_ms;
    out.observation_design_rank=design_rank;
    std::vector<double> per_atom_signal_squares(scenario.halos.size());
    std::vector<std::size_t> per_atom_signal_counts(scenario.halos.size());
    for (std::size_t h=0;h<scenario.halos.size();++h)
    {
        const auto & halo=scenario.halos[h];
        const auto & estimate=estimates[h];
        if (estimate.fallback) ++out.fallback_count;
        out.donor_total+=estimate.donor_count;
        if (!estimate.point) continue;
        if (!std::isfinite(estimate.point->GetAmplitude()) || estimate.point->GetAmplitude()<0 ||
            !std::isfinite(estimate.point->GetWidth()) || estimate.point->GetWidth()<=0 ||
            !std::isfinite(estimate.point->GetOffset()))
        {++out.invalid_estimates; continue;}
        ++out.available;
        const auto a_error=estimate.point->GetAmplitude()-halo.truth.GetAmplitude();
        const auto b_error=estimate.point->GetWidth()-halo.truth.GetWidth();
        const auto c_error=estimate.point->GetOffset()-halo.truth.GetOffset();
        out.parameter_a_square_sum+=a_error*a_error;
        out.parameter_b_square_sum+=b_error*b_error;
        out.parameter_c_square_sum+=c_error*c_error;
    }
    for (std::size_t node=0;node<scenario.raw_map_positions.size();++node)
    {
        const auto & position=scenario.raw_map_positions[node];
        const double observed_halo_signal=scenario.raw_map_values[node]-FixedTargetSignal(scenario,position);
        double predicted_halo_signal{};
        bool available=true;
        for (std::size_t h=0;h<scenario.halos.size();++h)
        {
            if (SquaredDistance(position,scenario.halos[h].position)>6.25) continue;
            const double true_signal=Contribution(
                scenario.halos[h].truth,scenario.halos[h].position,position);
            if (!estimates[h].point) {available=false; continue;}
            const double predicted=Contribution(*estimates[h].point,scenario.halos[h].position,position);
            predicted_halo_signal+=predicted;
            const auto per_atom_error=predicted-true_signal;
            per_atom_signal_squares[h]+=per_atom_error*per_atom_error;
            ++per_atom_signal_counts[h];
        }
        if (!available) {++out.unavailable_nodes; continue;}
        const auto error=predicted_halo_signal-observed_halo_signal;
        out.signal_square_sum+=error*error;
        ++out.signal_count;
    }
    std::vector<double> per_atom_rmses;
    for (std::size_t h=0;h<scenario.halos.size();++h)
        if (per_atom_signal_counts[h])
            per_atom_rmses.push_back(std::sqrt(per_atom_signal_squares[h]/
                static_cast<double>(per_atom_signal_counts[h])));
    if (!per_atom_rmses.empty())
    {
        std::sort(per_atom_rmses.begin(),per_atom_rmses.end());
        const auto middle=per_atom_rmses.size()/2;
        out.per_atom_signal_rmse_median=per_atom_rmses.size()%2 ? per_atom_rmses[middle] :
            (per_atom_rmses[middle-1]+per_atom_rmses[middle])/2.0;
        out.per_atom_signal_rmse_max=per_atom_rmses.back();
    }
    return out;
}

Summary Evaluate(const Policy & policy, Scenario & scenario)
{
    const auto start=std::chrono::steady_clock::now();
    auto tree=KDTreeAlgorithm<Donor>::BuildKDTree(scenario.donors);
    std::vector<Estimate> estimates;
    estimates.reserve(scenario.halos.size());
    for (const auto & halo:scenario.halos)
        estimates.push_back(EstimateHalo(policy,halo,scenario.donors,tree.get()));
    const auto elapsed=std::chrono::duration<double,std::milli>(
        std::chrono::steady_clock::now()-start).count();
    return ScoreEstimates(scenario,estimates,elapsed);
}

ObservationFit FitObservationDriven(Scenario & scenario)
{
    using Matrix=Eigen::MatrixXd;
    using Vector=Eigen::VectorXd;
    const Policy width_policy{"knn-4-unbounded",false,4,0.0,false};
    auto tree=KDTreeAlgorithm<Donor>::BuildKDTree(scenario.donors);
    ObservationFit out;
    out.estimates.reserve(scenario.halos.size());
    for (const auto & halo:scenario.halos)
        out.estimates.push_back(EstimateHalo(width_policy,halo,scenario.donors,tree.get()));

    const auto rows=static_cast<Eigen::Index>(scenario.raw_map_positions.size());
    const auto cols=static_cast<Eigen::Index>(2*scenario.halos.size());
    Matrix design(rows,cols);
    Vector response(rows);
    for (Eigen::Index row=0;row<rows;++row)
    {
        const auto & position=scenario.raw_map_positions[static_cast<std::size_t>(row)];
        response(row)=scenario.raw_map_values[static_cast<std::size_t>(row)]-
            FixedTargetSignal(scenario,position);
        for (std::size_t h=0;h<scenario.halos.size();++h)
        {
            const auto & halo=scenario.halos[h];
            const auto & point=*out.estimates[h].point;
            const auto basis=joint::EvaluateKernel(
                SquaredDistance(position,halo.position),point.GetWidth(),2.5);
            design(row,static_cast<Eigen::Index>(2*h))=basis.gaussian;
            design(row,static_cast<Eigen::Index>(2*h+1))=basis.charge;
        }
    }
    Eigen::ColPivHouseholderQR<Matrix> full_qr(design);
    full_qr.setThreshold(1e-10);
    out.design_rank=static_cast<std::size_t>(full_qr.rank());
    if (out.design_rank<static_cast<std::size_t>(cols))
    {
        for (auto & estimate:out.estimates) estimate.point.reset();
        return out;
    }
    double best_rss=std::numeric_limits<double>::infinity();
    std::vector<double> best_amplitudes(scenario.halos.size());
    std::vector<double> best_charges(scenario.halos.size());
    bool found_feasible{};
    const std::size_t active_sets=std::size_t{1}<<scenario.halos.size();
    for (std::size_t active=0;active<active_sets;++active)
    {
        std::vector<Eigen::Index> selected_columns;
        for (std::size_t h=0;h<scenario.halos.size();++h)
        {
            if ((active&(std::size_t{1}<<h))==0)
                selected_columns.push_back(static_cast<Eigen::Index>(2*h));
            selected_columns.push_back(static_cast<Eigen::Index>(2*h+1));
        }
        Matrix reduced(rows,static_cast<Eigen::Index>(selected_columns.size()));
        for (std::size_t col=0;col<selected_columns.size();++col)
            reduced.col(static_cast<Eigen::Index>(col))=design.col(selected_columns[col]);
        Eigen::ColPivHouseholderQR<Matrix> qr(reduced);
        qr.setThreshold(1e-10);
        if (qr.rank()!=reduced.cols()) continue;
        const Vector beta=qr.solve(response);
        if (!beta.allFinite()) continue;
        std::vector<double> amplitudes(scenario.halos.size());
        std::vector<double> charges(scenario.halos.size());
        Eigen::Index column{};
        bool feasible=true;
        for (std::size_t h=0;h<scenario.halos.size();++h)
        {
            if ((active&(std::size_t{1}<<h))==0)
            {
                amplitudes[h]=beta(column++);
                if (amplitudes[h]<0.0) feasible=false;
            }
            charges[h]=beta(column++);
        }
        if (!feasible) continue;
        const double rss=(reduced*beta-response).squaredNorm();
        if (rss<best_rss)
        {
            best_rss=rss;
            best_amplitudes=std::move(amplitudes);
            best_charges=std::move(charges);
            found_feasible=true;
        }
    }
    if (found_feasible)
        for (std::size_t h=0;h<scenario.halos.size();++h)
        {
            const auto width=out.estimates[h].point->GetWidth();
            out.estimates[h].point=GaussianModel3D{best_amplitudes[h],width,best_charges[h]};
        }
    else
        for (auto & estimate:out.estimates) estimate.point.reset();
    return out;
}

Summary EvaluateObservationDriven(Scenario & scenario)
{
    const auto start=std::chrono::steady_clock::now();
    const auto fit=FitObservationDriven(scenario);
    const auto elapsed=std::chrono::duration<double,std::milli>(
        std::chrono::steady_clock::now()-start).count();
    return ScoreEstimates(scenario,fit.estimates,elapsed,fit.design_rank);
}

template<class Evaluation>
Summary MedianTimedEvaluation(const Scenario & scenario, Evaluation evaluation)
{
    auto warmup=scenario;
    (void)evaluation(warmup);
    std::vector<double> timings;
    timings.reserve(7);
    Summary result;
    for (int repeat=0;repeat<7;++repeat)
    {
        auto input=scenario;
        result=evaluation(input);
        timings.push_back(result.cpu_ms);
    }
    std::sort(timings.begin(),timings.end());
    result.cpu_ms=timings[timings.size()/2];
    return result;
}

bool SameObservationDrivenAfterPermutation(const Scenario & scenario)
{
    auto forward=scenario;
    auto reverse=scenario;
    std::reverse(reverse.donors.begin(),reverse.donors.end());
    const auto a=FitObservationDriven(forward);
    const auto b=FitObservationDriven(reverse);
    if (a.design_rank!=b.design_rank || a.estimates.size()!=b.estimates.size()) return false;
    for (std::size_t h=0;h<a.estimates.size();++h)
    {
        if (a.estimates[h].point.has_value()!=b.estimates[h].point.has_value()) return false;
        if (a.estimates[h].point && (a.estimates[h].point->GetAmplitude()!=b.estimates[h].point->GetAmplitude() ||
            a.estimates[h].point->GetWidth()!=b.estimates[h].point->GetWidth() ||
            a.estimates[h].point->GetOffset()!=b.estimates[h].point->GetOffset())) return false;
    }
    return true;
}

bool SameEstimatesAfterPermutation(const Policy & policy, const Scenario & scenario)
{
    auto forward=scenario;
    auto reverse=scenario;
    std::reverse(reverse.donors.begin(),reverse.donors.end());
    auto forward_tree=KDTreeAlgorithm<Donor>::BuildKDTree(forward.donors);
    auto reverse_tree=KDTreeAlgorithm<Donor>::BuildKDTree(reverse.donors);
    for (std::size_t i=0;i<scenario.halos.size();++i)
    {
        const auto a=EstimateHalo(policy,scenario.halos[i],forward.donors,forward_tree.get());
        const auto b=EstimateHalo(policy,scenario.halos[i],reverse.donors,reverse_tree.get());
        if (a.donor_count!=b.donor_count || a.fallback!=b.fallback ||
            a.point.has_value()!=b.point.has_value()) return false;
        if (a.point && (a.point->GetAmplitude()!=b.point->GetAmplitude() ||
            a.point->GetWidth()!=b.point->GetWidth() || a.point->GetOffset()!=b.point->GetOffset())) return false;
    }
    return true;
}

std::optional<double> Rmse(double square_sum, std::size_t count)
{
    if (!count) return std::nullopt;
    return std::sqrt(square_sum/static_cast<double>(count));
}

void PrintMetric(std::optional<double> value)
{
    if (!value) std::cout << "NA";
    else std::cout << std::scientific << std::setprecision(6) << *value;
}

std::size_t WorkspaceBytes(const Scenario & scenario, bool observation_driven=false)
{
    const auto observation_workspace=observation_driven ? scenario.raw_map_positions.size()*
        (2*scenario.halos.size()*sizeof(double)+sizeof(double)) : 0;
    return scenario.donors.size()*(sizeof(Donor)+sizeof(Donor *)+sizeof(KDNode<Donor>))+
        scenario.halos.size()*(sizeof(Halo)+sizeof(Estimate))+
        scenario.donors.size()*sizeof(Donor *)+observation_workspace;
}
}

int main()
{
    std::cout << "scenario,policy,available_halos,total_halos,unavailable_nodes,invalid_estimates,"
        "observation_design_rank,mean_donors,fallback_halos,halo_signal_rmse,"
        "per_atom_signal_rmse_median,per_atom_signal_rmse_max,A_rmse,B_rmse,C_rmse,"
        "cpu_ms,deterministic,workspace_bytes\n";
    for (std::size_t scenario_index=0;scenario_index<10;++scenario_index)
    {
        for (const auto & policy:Policies())
        {
            auto scenario=MakeScenario(scenario_index);
            const auto result=MedianTimedEvaluation(scenario,
                [&policy](Scenario & input){return Evaluate(policy,input);});
            const auto deterministic=SameEstimatesAfterPermutation(policy,scenario);
            if (!deterministic) return 2;
            std::cout << scenario.name << ',' << policy.name << ',' << result.available << ',' << result.halos << ','
                << result.unavailable_nodes << ',' << result.invalid_estimates << ','
                << result.observation_design_rank << ','
                << std::fixed << std::setprecision(3)
                << static_cast<double>(result.donor_total)/static_cast<double>(result.halos) << ','
                << result.fallback_count << ',';
            PrintMetric(Rmse(result.signal_square_sum,result.signal_count));
            std::cout << ',';
            PrintMetric(result.per_atom_signal_rmse_median);
            std::cout << ',';
            PrintMetric(result.per_atom_signal_rmse_max);
            std::cout << ',';
            PrintMetric(Rmse(result.parameter_a_square_sum,result.available));
            std::cout << ',';
            PrintMetric(Rmse(result.parameter_b_square_sum,result.available));
            std::cout << ',';
            PrintMetric(Rmse(result.parameter_c_square_sum,result.available));
            std::cout << ',' << std::scientific << std::setprecision(6) << result.cpu_ms << ','
                << (deterministic ? "true" : "false") << ',' << WorkspaceBytes(scenario) << '\n';
        }
        auto scenario=MakeScenario(scenario_index);
        const auto result=MedianTimedEvaluation(scenario,
            [](Scenario & input){return EvaluateObservationDriven(input);});
        const auto deterministic=SameObservationDrivenAfterPermutation(scenario);
        if (!deterministic) return 3;
        std::cout << scenario.name << ",observation-ac-knn4-width,"
            << result.available << ',' << result.halos << ',' << result.unavailable_nodes << ','
            << result.invalid_estimates << ',' << result.observation_design_rank << ','
            << std::fixed << std::setprecision(3)
            << static_cast<double>(result.donor_total)/static_cast<double>(result.halos) << ','
            << result.fallback_count << ',';
        PrintMetric(Rmse(result.signal_square_sum,result.signal_count));
        std::cout << ',';
        PrintMetric(result.per_atom_signal_rmse_median);
        std::cout << ',';
        PrintMetric(result.per_atom_signal_rmse_max);
        std::cout << ',';
        PrintMetric(Rmse(result.parameter_a_square_sum,result.available));
        std::cout << ',';
        PrintMetric(Rmse(result.parameter_b_square_sum,result.available));
        std::cout << ',';
        PrintMetric(Rmse(result.parameter_c_square_sum,result.available));
        std::cout << ',' << std::scientific << std::setprecision(6) << result.cpu_ms << ','
            << (deterministic ? "true" : "false") << ',' << WorkspaceBytes(scenario,true) << '\n';
    }
}
