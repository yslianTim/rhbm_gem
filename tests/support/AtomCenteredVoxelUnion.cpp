#include "support/AtomCenteredVoxelUnion.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

namespace second_stage_test::matched::atom_union {
namespace {
template<class Visit>
void Sphere(const Atom & atom,const rhbm_gem::MapObject & map,double radius,Visit visit)
{
    if (!(radius>0) || !std::isfinite(radius)) throw std::invalid_argument("Invalid voxel radius.");
    const auto n=map.GetGridSize(); const auto origin=map.GetOrigin(),spacing=map.GetGridSpacing();
    std::array<int,3> lo{},hi{};
    for (std::size_t k=0;k<3;++k)
    {
        lo[k]=static_cast<int>(std::clamp(std::floor((atom.position[k]-radius-origin[k])/spacing[k]),0.0,static_cast<double>(n[k])));
        hi[k]=static_cast<int>(std::clamp(std::floor((atom.position[k]+radius-origin[k])/spacing[k]),-1.0,static_cast<double>(n[k]-1)));
    }
    for (int z=lo[2];z<=hi[2];++z) for (int y=lo[1];y<=hi[1];++y) for (int x=lo[0];x<=hi[0];++x)
    {
        const auto i=static_cast<std::size_t>(x)+static_cast<std::size_t>(n[0])*(static_cast<std::size_t>(y)+static_cast<std::size_t>(n[1])*static_cast<std::size_t>(z));
        const auto p=map.GetGridPosition(i); const double square=SquareDistance(p,atom.position);
        if (square<=radius*radius) visit(i,p,square);
    }
}
}
unique_grid::Grid BuildGrid(const std::vector<Atom> & atoms,const std::vector<Stencil> & stencils,
    const rhbm_gem::MapObject & generation,const rhbm_gem::MapObject & observation,double cutoff)
{
    if (generation.GetGridSize()!=observation.GetGridSize()) throw std::invalid_argument("Grid dimensions differ.");
    const auto size=generation.GetMapValueArraySize();
    std::vector<std::size_t> coverage(size); std::vector<double> nearest(size,std::numeric_limits<double>::infinity());
    for (const auto & atom:atoms) Sphere(atom,generation,cutoff,[&](std::size_t i,const Position &,double square)
    {++coverage[i]; nearest[i]=std::min(nearest[i],square);});
    unique_grid::Grid out; std::vector<int> rows(size,-1);
    out.voxels.reserve(static_cast<std::size_t>(std::count_if(coverage.begin(),coverage.end(),[](auto n){return n>0;})));
    for (std::size_t i=0;i<size;++i) if (coverage[i])
    {
        const double y=observation.GetMapValue(i); if (!std::isfinite(y)) throw std::invalid_argument("Nonfinite voxel observation.");
        rows[i]=static_cast<int>(out.voxels.size());
        out.voxels.push_back({i,coverage[i],generation.GetGridPosition(i),y,std::sqrt(nearest[i]),false});
    }
    for (const auto & stencil:stencils)
    {
        auto & slots=out.sample_rows.emplace_back();
        for (std::size_t k=0;k<64;++k)
        {
            const auto i=stencil.slots[k].index;
            if (i>=size || rows[i]<0) throw std::invalid_argument("Stencil voxel outside atom union.");
            slots[k]=static_cast<std::size_t>(rows[i]); out.voxels[slots[k]].in_stencil=true;
        }
    }
    return out;
}
joint_ac::Blocks BuildBlocks(const unique_grid::Grid & grid,const std::vector<Atom> & atoms,
    const rhbm_gem::MapObject & generation,const std::vector<double> & alphas,double cutoff)
{
    if (alphas.size()!=atoms.size()) throw std::invalid_argument("Block alpha population differs.");
    std::vector<int> rows(generation.GetMapValueArraySize(),-1);
    for (std::size_t p=0;p<grid.voxels.size();++p) rows.at(grid.voxels[p].index)=static_cast<int>(p);
    joint_ac::Blocks blocks;
    for (std::size_t a=0;a<atoms.size();++a)
    {
        joint_ac::Block block{a,alphas[a],{}};
        Sphere(atoms[a],generation,cutoff,[&](std::size_t i,const Position &,double)
        {
            if (rows[i]<0) throw std::invalid_argument("Block voxel outside union.");
            block.rows.push_back(rows[i]);
        });
        blocks.push_back(std::move(block));
    }
    return blocks;
}
Eigen::SparseMatrix<double> BuildDesign(const unique_grid::Grid & grid,const std::vector<Atom> & atoms,
    const rhbm_gem::MapObject & generation,double cutoff)
{
    std::vector<int> rows(generation.GetMapValueArraySize(),-1); std::size_t memberships{};
    for (std::size_t p=0;p<grid.voxels.size();++p) {rows.at(grid.voxels[p].index)=static_cast<int>(p); memberships+=grid.voxels[p].multiplicity;}
    Eigen::SparseMatrix<double> x(static_cast<Eigen::Index>(grid.voxels.size()),static_cast<Eigen::Index>(2*atoms.size()));
    x.reserve(static_cast<Eigen::Index>(2*memberships));
    for (std::size_t a=0;a<atoms.size();++a) for (int kind=0;kind<2;++kind)
    {
        x.startVec(static_cast<Eigen::Index>(2*a+static_cast<std::size_t>(kind)));
        Sphere(atoms[a],generation,cutoff,[&](std::size_t i,const Position &,double square)
        {
            if (rows[i]<0) return;
            const auto basis=EvaluateBasis(square,atoms[a].width,cutoff);
            const double value=kind==0 ? basis.gaussian : basis.charge;
            if (value!=0) x.insertBack(rows[i],static_cast<Eigen::Index>(2*a+static_cast<std::size_t>(kind)))=value;
        });
    }
    x.finalize(); return x;
}
} // namespace second_stage_test::matched::atom_union
