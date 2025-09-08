#pragma once

#include <vector>
#include <array>
#include <tuple>
#include "DataObjectVisitorBase.hpp"

class SphereSampler;

class MapInterpolationVisitor : public DataObjectVisitorBase
{
    SphereSampler * m_sampler;
    std::array<float, 3> m_position, m_axis_vector;
    std::vector<std::tuple<float, float>> m_sampling_data_list;
    std::vector<std::tuple<float, std::array<float, 3>>> m_point_list;

public:
    explicit MapInterpolationVisitor(SphereSampler * sampler);
    ~MapInterpolationVisitor() = default;
    void VisitMapObject(MapObject * data_object) override;

    void SetPosition(const std::array<float, 3> & position) { m_position = position; };
    void SetAxisVector(const std::array<float, 3> & axis_vector) { m_axis_vector = axis_vector; };
    const std::vector<std::tuple<float, float>> & GetSamplingDataList(void) const { return m_sampling_data_list; }
    std::vector<std::tuple<float, float>> && TakeSamplingDataList(void);

private:
    float MakeInterpolationInMapObject(MapObject * data_object, const std::array<float, 3> & position);

};
