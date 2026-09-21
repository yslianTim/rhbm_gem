#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include "core/command/detail/SimulationGeometry.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include <memory>
#include <string>

namespace joint_partial_test {
namespace core=rhbm_gem::core;
struct Fixture
{
    std::unique_ptr<rhbm_gem::ModelObject> model;
    std::unique_ptr<rhbm_gem::MapObject> map;
    std::vector<double> a,b,c;
};
inline Fixture Make(const std::string & name)
{
    Fixture f;
    std::vector<double> x=name=="bridge" ? std::vector<double>{-2.8,2.8,0} :
        name=="weak" ? std::vector<double>{0,4.8} : std::vector<double>{0,1.2};
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for(std::size_t i=0;i<x.size();++i)
    {
        auto atom=std::make_unique<rhbm_gem::AtomObject>();
        atom->SetSerialID(static_cast<int>(i+1)); atom->SetElement(Element::CARBON);
        atom->SetPosition(x[i],0,0); atom->SetChainID("A"); atom->SetComponentID("ALA");
        atom->SetAtomID(i==0 ? "CA" : "CB");
        atoms.push_back(std::move(atom));
        const auto index=static_cast<double>(i);
        f.a.push_back(2+.3*index); f.b.push_back(name=="weak" && i==1 ? .3 : .5+.05*index); f.c.push_back(.2-.05*index);
    }
    f.model=std::make_unique<rhbm_gem::ModelObject>(std::move(atoms));
    f.model->SelectAtoms([&](const auto & atom) {return name=="all" || atom.GetSerialID()==1 || (name=="bridge" && atom.GetSerialID()==2);});
    const std::array<int,3> dims{41,25,25};
    const std::array<double,3> spacing{.3,.3,.3},origin{-6,-3.6,-3.6};
    auto values=std::make_unique<double[]>(41*25*25);
    for(int z=0;z<dims[2];++z) for(int y=0;y<dims[1];++y) for(int k=0;k<dims[0];++k)
    {
        const auto position=core::simulation::GridPosition({k,y,z},spacing,origin);
        double value=0;
        for(std::size_t i=0;i<x.size();++i)
        {
            const double square=core::simulation::SupportSquare(position,{x[i],0,0});
            const auto kernel=core::joint_component::EvaluateKernel(square,f.b[i],2.5);
            value+=f.a[i]*kernel.gaussian+f.c[i]*kernel.charge;
        }
        values[static_cast<std::size_t>(k+dims[0]*(y+dims[1]*z))]=value;
    }
    f.map=std::make_unique<rhbm_gem::MapObject>(dims,spacing,origin,std::move(values));
    return f;
}
}
