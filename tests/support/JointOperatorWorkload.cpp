#include "JointOperatorWorkload.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <bit>
#include <boost/hash2/sha2.hpp>

namespace second_stage_test {
std::string OperatorWorkloadHash(const rhbm_gem::core::JointProblemInput & input)
{
    // Versioned canonical stream: uint64 little endian lengths/indices, exact
    // IEEE double bits and length-prefixed identity bytes. No snapshot copy.
    boost::hash2::sha2_256 hash;
    const auto integer=[&](std::uint64_t value) {
        std::array<unsigned char,8> bytes{};
        for(unsigned k=0;k<8;++k) bytes[k]=static_cast<unsigned char>((value>>(8*k))&255);
        hash.update(bytes.data(),bytes.size());
    };
    const auto string=[&](const std::string & value) {integer(value.size()); hash.update(value.data(),value.size());};
    string("frozen-lattice-input-v1"); integer(input.atom_ids.size()); integer(input.row_ids.size());
    for(const auto & id:input.atom_ids) string(id);
    for(std::size_t row=0;row<input.row_ids.size();++row)
    {string(input.row_ids[row]); integer(std::bit_cast<std::uint64_t>(input.observations[row]));}
    for(const auto & support:input.support)
    {
        integer(support.size());
        for(const auto & v:support) {integer(v.row); integer(std::bit_cast<std::uint64_t>(v.squared_distance));}
    }
    return boost::hash2::to_string(hash.result());
}
rhbm_gem::core::JointProblemInput OperatorWorkload(const std::string & topology,int count)
{
    if((topology!="chain" && topology!="cube") || count<=0 || count>10000)
        throw std::invalid_argument("Expected chain/cube and 1..10000 atoms");
    using Voxel=std::array<int,3>; // z,y,x; lexicographic ordering.
    std::vector<Voxel> offsets,centers,rows;
    for(int z=-5;z<=5;++z) for(int y=-5;y<=5;++y) for(int x=-5;x<=5;++x)
        if(x*x+y*y+z*z<=25) offsets.push_back({z,y,x});
    int side=1; while(side*side*side<count) ++side;
    rows.reserve(static_cast<std::size_t>(count)*offsets.size());
    for(int a=0;a<count;++a)
    {
        const Voxel center=topology=="chain" ? Voxel{0,0,7*a} : Voxel{7*(a/(side*side)),7*((a/side)%side),7*(a%side)};
        centers.push_back(center);
        for(const auto & d:offsets) rows.push_back({center[0]+d[0],center[1]+d[1],center[2]+d[2]});
    }
    std::sort(rows.begin(),rows.end()); rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
    rhbm_gem::core::JointProblemInput input;
    input.support.resize(static_cast<std::size_t>(count)); input.observations.resize(rows.size());
    input.row_ids.reserve(rows.size());
    for(std::size_t r=0;r<rows.size();++r)
    {
        const auto & v=rows[r]; input.row_ids.push_back(std::to_string(v[0])+":"+std::to_string(v[1])+":"+std::to_string(v[2]));
        input.observations[r]=1e-5*std::sin(.13*v[0]+.17*v[1]+.19*v[2]);
    }
    for(int a=0;a<count;++a)
    {
        input.atom_ids.push_back(std::to_string(a+1)); auto & support=input.support[static_cast<std::size_t>(a)];
        support.reserve(offsets.size()); const auto & center=centers[static_cast<std::size_t>(a)];
        for(const auto & d:offsets)
        {
            const Voxel v{center[0]+d[0],center[1]+d[1],center[2]+d[2]};
            const auto row=static_cast<std::size_t>(std::lower_bound(rows.begin(),rows.end(),v)-rows.begin());
            const double square=.25*(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
            support.push_back({row,square});
            const auto b=rhbm_gem::core::joint_component::EvaluateKernel(square,.5,2.5);
            input.observations[row]+=2*b.gaussian+.2*b.charge;
        }
    }
    return input;
}
}
