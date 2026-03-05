#pragma once

#include <vector>
#include <array>
#include <tuple>
#include "DataObjectVisitor.hpp"

class SamplerBase;

namespace rhbm_gem {

class MapInterpolationVisitor : public ConstDataObjectVisitor
{
    ::SamplerBase * m_sampler;
    std::array<float, 3> m_position, m_axis_vector;
    std::vector<std::tuple<float, float>> m_sampling_data_list;

public:
    explicit MapInterpolationVisitor(::SamplerBase * sampler);
    ~MapInterpolationVisitor() = default;
    void VisitAtomObject(const AtomObject & data_object) override;
    void VisitBondObject(const BondObject & data_object) override;
    void VisitModelObject(const ModelObject & data_object) override;
    void VisitMapObject(const MapObject & data_object) override;

    void SetPosition(const std::array<float, 3> & position) { m_position = position; };
    void SetAxisVector(const std::array<float, 3> & axis_vector) { m_axis_vector = axis_vector; };
    const std::vector<std::tuple<float, float>> & GetSamplingDataList() const { return m_sampling_data_list; }
    std::vector<std::tuple<float, float>> ConsumeSamplingDataList();

private:
    float MakeInterpolationInMapObject(
        const MapObject & data_object, const std::array<float, 3> & position);

};

} // namespace rhbm_gem
