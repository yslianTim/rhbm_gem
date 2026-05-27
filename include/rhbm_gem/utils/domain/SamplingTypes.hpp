#pragma once

#include <array>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <vector>

enum class SphereSamplingMethod
{
    RadiusUniformRandom,
    VolumeUniformRandom,
    FibonacciDeterministic
};

struct SamplingPoint
{
    float distance{ 0.0f };
    std::array<float, 3> position{ 0.0f, 0.0f, 0.0f };
    bool is_selected{ true };
};

using SamplingPointList = std::vector<SamplingPoint>;

struct LocalPotentialSample
{
    float response{ 0.0f };
    SamplingPoint point{};
};

using LocalPotentialSampleList = std::vector<LocalPotentialSample>;

struct SeriesPoint
{
    std::vector<double> basis_list{};
    double response{ 0.0 };

    SeriesPoint() = default;

    SeriesPoint(
        std::initializer_list<double> basis_list_value,
        double response_value) :
        basis_list{ basis_list_value },
        response{ response_value }
    {
    }

    SeriesPoint(
        std::vector<double> basis_list_value,
        double response_value) :
        basis_list{ std::move(basis_list_value) },
        response{ response_value }
    {
    }

    size_t GetBasisSize() const noexcept { return basis_list.size(); }

    double GetBasisValue(size_t index) const
    {
        if (index >= basis_list.size())
        {
            throw std::out_of_range("SeriesPoint basis value index is out of range.");
        }
        return basis_list.at(index);
    }
};

using SeriesPointList = std::vector<SeriesPoint>;
