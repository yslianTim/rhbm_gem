#pragma once
#include "core/detail/joint_component/Numerics.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <cmath>
#include <memory>
#include <vector>

namespace second_stage_test {
// Known normalized Gaussian truth, with deterministic geometry/noise phases.
inline std::unique_ptr<rhbm_gem::ModelObject> FirstStageSyntheticModel(double separation, int replicate)
{
    using namespace rhbm_gem;
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for (int i=0;i<24;++i)
    {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(i+1); atom->SetElement(Element::CARBON);
        atom->SetPosition(.13+.07*replicate+separation*(i%4),.17-.03*replicate+separation*((i/4)%3),
            .19+.11*replicate+separation*(i/12));
        atoms.push_back(std::move(atom));
    }
    auto model=std::make_unique<ModelObject>(std::move(atoms));
    model->SelectAllAtoms(); return model;
}
inline std::unique_ptr<rhbm_gem::MapObject> FirstStageSyntheticMap(const rhbm_gem::ModelObject & model,double spacing,double noise,int replicate)
{
    using namespace rhbm_gem;
    namespace n = rhbm_gem::core::joint_component;
    const auto last=model.GetAtomList().back()->GetPosition();
    const std::array<double,3> origin{-3.,-3.,-3.};
    std::array<int,3> dims{};
    for (std::size_t k=0;k<3;++k) dims[k]=static_cast<int>(std::ceil((last[k]+6.)/spacing));
    const auto count=static_cast<std::size_t>(dims[0])*static_cast<std::size_t>(dims[1])*static_cast<std::size_t>(dims[2]);
    auto values=std::make_unique<double[]>(count);
    for (int z=0;z<dims[2];++z) for (int y=0;y<dims[1];++y) for (int x=0;x<dims[0];++x)
    {
        const std::array<double,3> p{origin[0]+x*spacing,origin[1]+y*spacing,origin[2]+z*spacing};
        double value=noise*std::sin(2.13*x+1.57*y+.73*z+.37*replicate);
        for (std::size_t i=0;i<model.GetAtomList().size();++i)
        {
            const auto q=model.GetAtomList()[i]->GetPosition(); double square=0;
            for (std::size_t k=0;k<3;++k) square+=(p[k]-q[k])*(p[k]-q[k]);
            value+=(1.5+.02*static_cast<double>(i))*n::EvaluateKernel(square,.35+.005*static_cast<double>(i),2.5).gaussian;
        }
        const auto index=static_cast<std::size_t>(x)+static_cast<std::size_t>(dims[0])*(static_cast<std::size_t>(y)+static_cast<std::size_t>(dims[1])*static_cast<std::size_t>(z));
        values[index]=value;
    }
    return std::make_unique<MapObject>(dims,std::array<double,3>{spacing,spacing,spacing},origin,std::move(values));
}
} // namespace second_stage_test
