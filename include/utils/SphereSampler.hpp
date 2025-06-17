#pragma once

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <array>
#include <vector>
#include <tuple>

class SphereSampler
{
    unsigned int m_thread_size, m_sampling_size;
    double m_distance_min, m_distance_max;
    //std::vector<std::tuple<float, std::array<float, 3>>> m_sampling_position_list;

public:
    SphereSampler(void);
    ~SphereSampler();

    void Print(void) const;
    //const std::vector<std::tuple<float, std::array<float, 3>>> & GenerateSamplingPoints(std::array<float, 3> position);
    std::vector<std::tuple<float, std::array<float, 3>>> GenerateSamplingPoints(std::array<float, 3> position);
    void SetThreadSize(unsigned int value) { m_thread_size = value; }
    void SetSamplingSize(unsigned int value) { m_sampling_size = value; }
    void SetDistanceRangeMinimum(double value) { m_distance_min = value; }
    void SetDistanceRangeMaximum(double value) { m_distance_max = value; }

private:


};