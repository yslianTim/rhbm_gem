#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
using Position=std::array<double,3>;
struct Halo
{
    int id{};
    Position position{};
    std::vector<std::size_t> rows;
    const rhbm_gem::AtomObject * atom{};
};
struct ClosePair
{
    const Halo * left{};
    const Halo * right{};
    double distance{};
};
class CoutCapture
{
public:
    CoutCapture():m_previous(std::cout.rdbuf(m_sink.rdbuf())) {}
    ~CoutCapture() {Restore();}
    void Restore()
    {
        if(m_previous) {std::cout.rdbuf(m_previous);m_previous=nullptr;}
    }
private:
    std::ostringstream m_sink;
    std::streambuf * m_previous{};
};
struct DistanceBin
{
    const char * name{};
    double lower{},upper{};
    std::size_t pairs{},overlap_pairs{},shared_voxels{};
    std::vector<std::size_t> shared_counts;
};

double Distance(const Position & left,const Position & right)
{
    double square{};
    for(std::size_t axis=0;axis<3;++axis)
    {
        const double delta=left[axis]-right[axis];
        square+=delta*delta;
    }
    return std::sqrt(square);
}

double Quantile(std::vector<double> values,double probability)
{
    if(values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(),values.end());
    const double index=probability*static_cast<double>(values.size()-1);
    const auto lower=static_cast<std::size_t>(std::floor(index));
    const auto upper=static_cast<std::size_t>(std::ceil(index));
    return values[lower]+(index-static_cast<double>(lower))*(values[upper]-values[lower]);
}

std::size_t SharedRows(const Halo & left,const Halo & right)
{
    std::size_t shared{};
    auto a=left.rows.begin(),b=right.rows.begin();
    while(a!=left.rows.end() && b!=right.rows.end())
    {
        if(*a<*b) ++a;
        else if(*b<*a) ++b;
        else {++shared;++a;++b;}
    }
    return shared;
}

void PrintDistanceSummary(const std::vector<double> & values)
{
    std::cout<<"{\"count\": "<<values.size();
    if(!values.empty())
        std::cout<<", \"minimum\": "<<*std::min_element(values.begin(),values.end())
            <<", \"p05\": "<<Quantile(values,0.05)<<", \"median\": "<<Quantile(values,0.5)
            <<", \"p95\": "<<Quantile(values,0.95)<<", \"maximum\": "
            <<*std::max_element(values.begin(),values.end());
    std::cout<<"}";
}
}

int main(int argc,char ** argv)
try
{
    if(argc!=4)
    {
        std::cerr<<"usage: joint_halo_geometry_census <map> <model.cif> <chain-id>\n";
        return 2;
    }
    CoutCapture parser_log_capture;
    auto map=rhbm_gem::ReadMap(argv[1]);
    auto model=rhbm_gem::ReadModel(argv[2]);
    if(!map || !model) throw std::runtime_error("Could not read the map or model.");
    const std::string chain=argv[3];
    model->SelectAtoms([&](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()==chain;
    });
    const auto problem=rhbm_gem::core::BuildPureHaloSupportExclusionProblem(*map,*model);
    const auto & input=problem.OriginalInput();
    if(!input.selection_domain) throw std::runtime_error("OriginalInput has no selection domain.");
    if(std::set<std::string>(input.row_ids.begin(),input.row_ids.end()).size()!=input.row_ids.size())
        throw std::runtime_error("OriginalInput row identities are not unique.");

    std::unordered_set<int> target_ids;
    for(const auto index:input.selection_domain->target_indices)
        target_ids.insert(std::stoi(input.atom_ids.at(index)));
    std::unordered_map<int,std::size_t> original_contributors;
    original_contributors.reserve(input.atom_ids.size());
    for(std::size_t index=0;index<input.atom_ids.size();++index)
        original_contributors.emplace(std::stoi(input.atom_ids[index]),index);

    std::vector<Halo> halos;
    std::unordered_set<int> halo_ids;
    for(const auto & atom:model->GetAtomList())
    {
        if(atom->GetElement()==Element::HYDROGEN || target_ids.contains(atom->GetSerialID())) continue;
        Halo halo{atom->GetSerialID(),atom->GetPosition(),{},atom.get()};
        if(const auto contributor=original_contributors.find(halo.id);contributor!=original_contributors.end())
        {
            for(const auto & support:input.support.at(contributor->second))
                halo.rows.push_back(support.row);
            std::sort(halo.rows.begin(),halo.rows.end());
            halo.rows.erase(std::unique(halo.rows.begin(),halo.rows.end()),halo.rows.end());
        }
        halo_ids.insert(halo.id);
        halos.push_back(std::move(halo));
    }
    std::sort(halos.begin(),halos.end(),[](const Halo & left,const Halo & right){return left.id<right.id;});
    std::vector<const rhbm_gem::AtomObject *> targets;
    targets.reserve(target_ids.size());
    for(const auto & atom:model->GetAtomList())
        if(target_ids.contains(atom->GetSerialID())) targets.push_back(atom.get());
    std::sort(targets.begin(),targets.end(),[](const auto * left,const auto * right){
        return left->GetSerialID()<right->GetSerialID();
    });

    std::map<std::size_t,std::size_t> observation_count_distribution;
    std::vector<double> halo_target_nearest_distances;
    std::vector<const Halo *> supported;
    std::size_t original_halo_affected_rows{};
    std::vector<bool> halo_affected(input.observations.size(),false);
    for(const auto & halo:halos)
    {
        ++observation_count_distribution[halo.rows.size()];
        if(!halo.rows.empty())
        {
            supported.push_back(&halo);
            for(const auto row:halo.rows) halo_affected.at(row)=true;
        }
        if(!targets.empty())
        {
            double nearest=std::numeric_limits<double>::infinity();
            for(const auto * target:targets) nearest=std::min(nearest,Distance(halo.position,target->GetPosition()));
            halo_target_nearest_distances.push_back(nearest);
        }
    }
    original_halo_affected_rows=static_cast<std::size_t>(std::count(halo_affected.begin(),halo_affected.end(),true));

    std::vector<double> supported_halo_nearest_distances;
    std::vector<double> all_halo_nearest_distances;
    std::vector<ClosePair> close_supported_pairs,close_halo_pairs;
    std::size_t supported_center_pairs_below_1_5{},all_center_pairs_below_1_5{},
        all_center_pairs_below_1_0{},all_center_pairs_below_0_5{},all_center_pairs_coincident{};
    for(const auto & halo:halos)
    {
        const auto * atom=model->FindAtomPtr(halo.id);
        const auto close=model->FindNeighborAtoms(*atom,1.5,false);
        for(const auto * neighbor:close)
        {
            if(neighbor->GetSerialID()<=halo.id || !halo_ids.contains(neighbor->GetSerialID())) continue;
            const double distance=Distance(halo.position,neighbor->GetPosition());
            if(distance<1.5)
            {
                ++all_center_pairs_below_1_5;
                const auto other=std::lower_bound(halos.begin(),halos.end(),neighbor->GetSerialID(),
                    [](const Halo & candidate,int id){return candidate.id<id;});
                close_halo_pairs.push_back({&halo,&*other,distance});
                if(distance<1.0) ++all_center_pairs_below_1_0;
                if(distance<0.5) ++all_center_pairs_below_0_5;
                if(distance<=1e-9) ++all_center_pairs_coincident;
            }
        }
        double radius=1.5,nearest=std::numeric_limits<double>::infinity();
        for(;;)
        {
            const auto neighbors=model->FindNeighborAtoms(*atom,radius,false);
            for(const auto * neighbor:neighbors)
                if(halo_ids.contains(neighbor->GetSerialID()))
                {
                    nearest=Distance(halo.position,neighbor->GetPosition());
                    break;
                }
            if(std::isfinite(nearest) || radius>=1024.0) break;
            radius*=2.0;
        }
        if(std::isfinite(nearest)) all_halo_nearest_distances.push_back(nearest);
    }
    std::array<DistanceBin,6> bins{{
        {"0_to_1_5_A",0.0,1.5,0,0,0,{}},{"1_5_to_2_A",1.5,2.0,0,0,0,{}},
        {"2_to_3_A",2.0,3.0,0,0,0,{}},{"3_to_4_A",3.0,4.0,0,0,0,{}},
        {"4_to_5_A",4.0,5.0,0,0,0,{}},
        {"5_A_or_more",5.0,std::numeric_limits<double>::infinity(),0,0,0,{}}}};
    for(std::size_t i=0;i<supported.size();++i)
    {
        double nearest=std::numeric_limits<double>::infinity();
        for(std::size_t j=0;j<supported.size();++j)
        {
            if(i==j) continue;
            const double distance=Distance(supported[i]->position,supported[j]->position);
            nearest=std::min(nearest,distance);
            if(j<i) continue;
            const auto shared=SharedRows(*supported[i],*supported[j]);
            for(auto & bin:bins) if(distance>=bin.lower && distance<bin.upper)
            {
                ++bin.pairs;
                if(shared)
                {
                    ++bin.overlap_pairs;
                    bin.shared_voxels+=shared;
                    bin.shared_counts.push_back(shared);
                }
                break;
            }
            if(distance<1.5)
            {
                ++supported_center_pairs_below_1_5;
                close_supported_pairs.push_back({supported[i],supported[j],distance});
            }
        }
        if(std::isfinite(nearest)) supported_halo_nearest_distances.push_back(nearest);
    }

    std::size_t fixed_count{},free_count{},support_overlap_pairs{};
    std::vector<double> shared_voxel_counts;
    for(const auto & bin:bins)
    {
        support_overlap_pairs+=bin.overlap_pairs;
        shared_voxel_counts.insert(shared_voxel_counts.end(),bin.shared_counts.begin(),bin.shared_counts.end());
    }
    for(const auto & [count,number]:observation_count_distribution)
        if(count<3) fixed_count+=number; else free_count+=number;
    const auto & exclusion=*input.selection_domain->pure_halo_support_exclusion;
    const auto by_distance=[](const ClosePair & left,const ClosePair & right) {
        if(left.distance!=right.distance) return left.distance<right.distance;
        if(left.left->id!=right.left->id) return left.left->id<right.left->id;
        return left.right->id<right.right->id;
    };
    std::sort(close_supported_pairs.begin(),close_supported_pairs.end(),by_distance);
    std::sort(close_halo_pairs.begin(),close_halo_pairs.end(),by_distance);
    parser_log_capture.Restore();
    std::cout<<std::setprecision(10)<<"{\n"
        <<"  \"dataset\": \"6Z6U chain "<<chain<<"\",\n"
        <<"  \"model\": \""<<argv[2]<<"\",\n"
        <<"  \"map\": \""<<argv[1]<<"\",\n"
        <<"  \"kernel_support_radius_angstrom\": 2.5,\n"
        <<"  \"target_count\": "<<target_ids.size()<<",\n"
        <<"  \"halo_candidates_in_model\": "<<halos.size()<<",\n"
        <<"  \"supported_halos\": "<<supported.size()<<",\n"
        <<"  \"original_observations\": "<<input.observations.size()<<",\n"
        <<"  \"original_row_ids_unique\": true,\n"
        <<"  \"original_row_identity\": \"OriginalInput row IDs map to unique global map voxel IDs\",\n"
        <<"  \"halo_affected_observations\": "<<original_halo_affected_rows<<",\n"
        <<"  \"zero_support_halos_from_problem_provenance\": "<<exclusion.zero_support_halos<<",\n"
        <<"  \"observation_count_distribution\": {";
    bool first=true;
    for(const auto & [count,number]:observation_count_distribution)
    {
        if(!first) std::cout<<", ";
        first=false;
        std::cout<<"\""<<count<<"\": "<<number;
    }
    std::cout<<"},\n  \"n_h_lt_3\": {\"count\": "<<fixed_count<<", \"fraction\": "
        <<(halos.empty()?0.0:static_cast<double>(fixed_count)/static_cast<double>(halos.size()))
        <<"},\n  \"n_h_ge_3\": {\"count\": "<<free_count<<", \"fraction\": "
        <<(halos.empty()?0.0:static_cast<double>(free_count)/static_cast<double>(halos.size()))<<"},\n"
        <<"  \"supported_halo_pairs\": "<<(supported.size()*(supported.size()-1)/2)<<",\n"
        <<"  \"all_halo_pairs_below_1_5_angstrom\": "<<all_center_pairs_below_1_5<<",\n"
        <<"  \"all_halo_pairs_below_1_0_angstrom\": "<<all_center_pairs_below_1_0<<",\n"
        <<"  \"all_halo_pairs_below_0_5_angstrom\": "<<all_center_pairs_below_0_5<<",\n"
        <<"  \"all_halo_coincident_center_pairs\": "<<all_center_pairs_coincident<<",\n"
        <<"  \"all_halo_nearest_center_distance_angstrom\": ";
    PrintDistanceSummary(all_halo_nearest_distances);
    std::cout<<",\n  \"supported_halo_pairs_below_1_5_angstrom\": "<<supported_center_pairs_below_1_5<<",\n"
        <<"  \"supported_halo_nearest_center_distance_angstrom\": ";
    PrintDistanceSummary(supported_halo_nearest_distances);
    std::cout<<",\n  \"halo_target_nearest_distance_angstrom\": ";
    PrintDistanceSummary(halo_target_nearest_distances);
    std::cout<<",\n  \"center_distance_bins_and_discrete_shared_voxels\": [\n";
    for(std::size_t index=0;index<bins.size();++index)
    {
        const auto & bin=bins[index];
        std::cout<<"    {\"distance_bin\": \""<<bin.name<<"\", \"pairs\": "<<bin.pairs
            <<", \"support_overlap_pairs\": "<<bin.overlap_pairs<<", \"shared_voxel_sum\": "<<bin.shared_voxels;
        if(!bin.shared_counts.empty())
        {
            std::vector<double> counts(bin.shared_counts.begin(),bin.shared_counts.end());
            std::cout<<", \"shared_voxels_per_overlapping_pair\": {\"median\": "<<Quantile(counts,0.5)
                <<", \"p95\": "<<Quantile(counts,0.95)<<", \"maximum\": "
                <<*std::max_element(bin.shared_counts.begin(),bin.shared_counts.end())<<"}";
        }
        std::cout<<"}"<<(index+1<bins.size()?",":"")<<"\n";
    }
    std::cout<<"  ],\n  \"support_overlap_pairs\": "<<support_overlap_pairs
        <<",\n  \"shared_voxels_per_overlap_pair\": {\"minimum\": "
        <<(shared_voxel_counts.empty()?0.0:*std::min_element(shared_voxel_counts.begin(),shared_voxel_counts.end()))
        <<", \"median\": "<<Quantile(shared_voxel_counts,0.5)<<", \"p95\": "
        <<Quantile(shared_voxel_counts,0.95)<<", \"maximum\": "
        <<(shared_voxel_counts.empty()?0.0:*std::max_element(shared_voxel_counts.begin(),shared_voxel_counts.end()))<<"}\n";
    const auto print_pairs=[&](const char * name,const std::vector<ClosePair> & pairs) {
        static bool first_pair_list=true;
        if(first_pair_list) std::cout<<",\n";
        else std::cout<<",\n";
        first_pair_list=false;
        std::cout<<"  \""<<name<<"\": [";
        const auto count=std::min<std::size_t>(pairs.size(),10);
        for(std::size_t index=0;index<count;++index)
        {
            const auto & pair=pairs[index];
            if(index) std::cout<<", ";
            std::cout<<"{\"left_atom_id\": "<<pair.left->id<<", \"right_atom_id\": "<<pair.right->id
                <<", \"distance_angstrom\": "<<pair.distance
                <<", \"left_atom_name\": \""<<pair.left->atom->GetAtomID()<<"\", \"right_atom_name\": \""
                <<pair.right->atom->GetAtomID()<<"\", \"left_chain\": \""<<pair.left->atom->GetChainID()
                <<"\", \"right_chain\": \""<<pair.right->atom->GetChainID()<<"\", \"left_component\": \""
                <<pair.left->atom->GetComponentID()<<"\", \"right_component\": \""<<pair.right->atom->GetComponentID()
                <<"\", \"left_element_atomic_number\": "<<static_cast<int>(pair.left->atom->GetElement())
                <<", \"right_element_atomic_number\": "<<static_cast<int>(pair.right->atom->GetElement())<<"}";
        }
        std::cout<<"]";
    };
    print_pairs("closest_supported_halo_pairs_below_1_5_angstrom",close_supported_pairs);
    print_pairs("closest_all_halo_pairs_below_1_5_angstrom",close_halo_pairs);
    std::cout<<"\n}\n";
    return 0;
}
catch(const std::exception & error)
{
    std::cerr<<"joint_halo_geometry_census: "<<error.what()<<'\n';
    return 1;
}
