#pragma once

#include <optional>
#include <string_view>

#include <Eigen/Dense>

namespace rhbm_gem
{

struct GaussianModel3DEvaluation
{
    double signal{ 0.0 };
    double offset_basis{ 0.0 };
    double response{ 0.0 };
};

class GaussianModel3D
{
    double m_amplitude{ 0.0 };
    double m_width{ 1.0 };
    double m_offset{ 0.0 };

    static constexpr int kParameterSize{ 3 };
    static constexpr int kAmplitudeIndex{ 0 };
    static constexpr int kWidthIndex{ 1 };
    static constexpr int kOffsetIndex{ 2 };
    static constexpr int kTransformedCoordinateSize{ 3 };
    static constexpr int kLogPeakHeightCoordinateIndex{ 0 };
    static constexpr int kLogWidthCoordinateIndex{ 1 };
    static constexpr int kOffsetToPeakRatioCoordinateIndex{ 2 };

public:
    using TransformedCoordinates = Eigen::Vector3d;

    GaussianModel3D() = default;
    GaussianModel3D(double amplitude, double width, double offset = 0.0);

    static constexpr int ParameterSize() { return kParameterSize; }
    static constexpr int AmplitudeIndex() { return kAmplitudeIndex; }
    static constexpr int WidthIndex() { return kWidthIndex; }
    static constexpr int OffsetIndex() { return kOffsetIndex; }
    static constexpr int TransformedCoordinateSize() { return kTransformedCoordinateSize; }
    static constexpr int LogPeakHeightCoordinateIndex() { return kLogPeakHeightCoordinateIndex; }
    static constexpr int LogWidthCoordinateIndex() { return kLogWidthCoordinateIndex; }
    static constexpr int OffsetToPeakRatioCoordinateIndex()
    {
        return kOffsetToPeakRatioCoordinateIndex;
    }

    static void RequireParameterVector(
        const Eigen::VectorXd & parameters,
        std::string_view value_name = "GaussianModel3D parameter vector");

    static GaussianModel3D FromVector(const Eigen::VectorXd & parameters);
    static GaussianModel3D FromVectorPrefix(const Eigen::VectorXd & parameters);
    static std::optional<GaussianModel3D> FromTransformedCoordinates(
        const TransformedCoordinates & coordinates);

    static void RequireFiniteModel(
        const GaussianModel3D & model,
        std::string_view value_name = "GaussianModel3D");
    static void RequireFinitePositiveWidthModel(
        const GaussianModel3D & model,
        std::string_view value_name = "GaussianModel3D");

    GaussianModel3D WithAmplitude(double value) const;
    GaussianModel3D WithWidth(double value) const;
    GaussianModel3D WithOffset(double value) const;

    double GetAmplitude() const { return m_amplitude; }
    double GetWidth() const { return m_width; }
    double GetOffset() const { return m_offset; }
    double GetHeight() const;
    Eigen::VectorXd ToVector() const;
    std::optional<TransformedCoordinates> ToTransformedCoordinates() const;
    double GetModelParameter(int par_id) const;
    double GetDisplayParameter(int par_id) const;
    double Intensity() const;
    GaussianModel3DEvaluation EvaluateAtDistance(double distance) const;
    double SignalAtDistance(double distance) const;
    double OffsetBasisAtDistance(double distance) const;
    double ResponseAtDistance(double distance) const;
};

class GaussianModel3DUncertainty
{
    double m_amplitude{ 0.0 };
    double m_width{ 0.0 };
    double m_offset{ 0.0 };

public:
    GaussianModel3DUncertainty() = default;
    GaussianModel3DUncertainty(double amplitude, double width, double offset = 0.0);

    static void RequireFiniteNonNegativeUncertainty(
        const GaussianModel3DUncertainty & uncertainty,
        std::string_view value_name = "GaussianModel3DUncertainty");

    double GetAmplitude() const { return m_amplitude; }
    double GetWidth() const { return m_width; }
    double GetOffset() const { return m_offset; }
    Eigen::VectorXd ToVector() const;
    double GetModelParameter(int par_id) const;
};

class GaussianModel3DWithUncertainty
{
    GaussianModel3D m_model{ 0.0, 0.0 };
    GaussianModel3DUncertainty m_standard_deviation{};

public:
    GaussianModel3DWithUncertainty() = default;
    GaussianModel3DWithUncertainty(
        GaussianModel3D model,
        GaussianModel3DUncertainty standard_deviation);

    const GaussianModel3D & GetModel() const { return m_model; }
    const GaussianModel3DUncertainty & GetStandardDeviationModel() const { return m_standard_deviation; }
    double GetDisplayParameter(int par_id) const;
    double GetDisplayStandardDeviation(int par_id) const;
    double GetModelParameter(int par_id) const;
    double GetModelStandardDeviation(int par_id) const;
    double IntensityStandardDeviation() const;
};

} // namespace rhbm_gem
