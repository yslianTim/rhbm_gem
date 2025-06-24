#pragma once

#include <vector>
#include <array>
#include <tuple>
#include "DataObjectVisitorBase.hpp"

class SphereSampler;

class MapInterpolationVisitor : public DataObjectVisitorBase
{
    SphereSampler * m_sphere_sampler;
    std::array<float, 3> m_position;
    std::vector<std::tuple<float, float>> m_sampling_data_list;
    std::vector<std::tuple<float, std::array<float, 3>>> m_sampling_point_list;

public:
    MapInterpolationVisitor(const std::vector<std::tuple<float, std::array<float, 3>>> & sampling_points);
    MapInterpolationVisitor(SphereSampler * sphere_sampler);
    ~MapInterpolationVisitor();
    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void Analysis(DataObjectManager * data_manager) override;

    void SetPosition(const std::array<float, 3> position) { m_position = position; };
    const std::vector<std::tuple<float, float>> & GetSamplingDataList(void) const { return m_sampling_data_list; }

private:
    float MakeInterpolationInMapObject(MapObject * data_object, std::array<float, 3> position);

};
