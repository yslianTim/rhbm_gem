#include "core/classifier/UmapEmbeddingOutput.hpp"
#include "core/classifier/UmapEmbeddingFeatures.hpp"
#include "core/classifier/UmapSpotCategories.hpp"

#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <boost/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

namespace rhbm_gem::core::classifier {

namespace {

bool WriteEmbedding(
    const std::filesystem::path & output_path,
    const PreparedUmapInput & prepared,
    const UmapEmbeddingResult & first,
    const UmapEmbeddingResult & second,
    std::string & error_message)
{
    std::ofstream output{ output_path, std::ios::out | std::ios::trunc };
    if (!output.is_open())
    {
        error_message = "Failed to open UMAP output file '"
            + output_path.string() + "'.";
        return false;
    }

    output << BuildLocalFittingCsvHeader()
        << ",umap x,umap y,mean shift x,mean shift y,bin number,second umap x,second umap y\n";
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t observation = 0;
        observation < prepared.feature_rows.size();
        ++observation)
    {
        const auto & row{ prepared.feature_rows[observation] };
        output << row.serial_id << ',' << row.residue << ',' << row.spot;
        for (std::size_t feature = 0;
            feature < kLocalFittingFeatureCount;
            ++feature)
        {
            output << ',';
            if (!std::isfinite(row.features[feature])) continue;
            output << row.features[feature];
        }
        const auto offset{ observation * kOutputDimensionCount };
        output << ',' << first.coordinates[offset] << ',' << first.coordinates[offset + 1]
            << ',' << first.mean_shift_coordinates[offset] << ',' << first.mean_shift_coordinates[offset + 1]
            << ',' << first.binning.point_bin_numbers[observation];
        if (second.coordinates.empty()) output << ",,";
        else output << ',' << second.coordinates[offset] << ',' << second.coordinates[offset + 1];
        output << '\n';
    }
    output.close();
    boost::json::array features;
    for (const auto name : GetUmapFeatureNames()) features.emplace_back(name);
    boost::json::array second_umap_features;
    for (const auto name : GetSecondUmapFeatureNames()) second_umap_features.emplace_back(name);
    boost::json::value map_normalization;
    if (const auto & normalization{ prepared.map_normalization }; normalization)
    {
        map_normalization = boost::json::object{
            { "requested", normalization->requested },
            { "applied", normalization->applied },
            { "divisor", normalization->divisor }
        };
    }
    boost::json::object exclusions;
    for (const auto & [reason, count] : prepared.exclusions) exclusions[reason] = count;

    std::ofstream metadata(output_path.string() + ".metadata.json");
    boost::json::object metadata_fields{
        { "estimator", "joint-components" },
        { "features", features },
        { "second_umap_features", second_umap_features },
        { "peeling_mode", "original-grid-v1" },
        { "map_normalization", map_normalization },
        { "column_semantics", boost::json::object{
            { kLocalFittingFeatureNames[4], "Joint Gaussian amplitude A" },
            { kLocalFittingFeatureNames[5], "Joint shared Gaussian/charge-kernel width tau" },
            { kLocalFittingFeatureNames[6], "Joint charge-kernel coefficient C" } } },
        { "standardization", "population-z-score; constant columns zero" },
        { "included_rows", prepared.feature_rows.size() },
        { "excluded_rows", prepared.excluded_rows },
        { "exclusions", exclusions }
    };
    metadata << boost::json::serialize(metadata_fields);
    metadata.close();
    if (!metadata)
    {
        error_message = "Failed to write UMAP metadata.";
        return false;
    }
    if (!output)
    {
        error_message = "Failed while writing UMAP output file '"
            + output_path.string() + "'.";
        return false;
    }
    return true;
}

bool WriteMeanShiftBins(
    const PreparedUmapInput & prepared,
    const UmapBinning & binning,
    const std::filesystem::path & output_path,
    bool include_element_ratios,
    std::string & error_message)
{
    std::ofstream output{ output_path };
    if (!output.is_open())
    {
        error_message = "Failed to open UMAP mean shift bin output file '" + output_path.string() + "'.";
        return false;
    }
    output << "bin number,bin x,bin y,point count,amplitude joint median,width joint median,offset joint median,"
        "spot class count,1st dominate spot,1st dominate spot count,1st dominate spot count ratio,"
        "2nd dominate spot,2nd dominate spot count,2nd dominate spot count ratio,"
        "3rd dominate spot,3rd dominate spot count,3rd dominate spot count ratio";
    if (include_element_ratios) output << ",carbon ratio,nitrogen ratio,oxygen ratio";
    output << '\n' << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto & [index, bin] : binning.bins)
    {
        std::map<std::string, std::size_t> spot_counts;
        std::size_t carbon_count{};
        std::size_t nitrogen_count{};
        std::size_t oxygen_count{};
        for (const auto point : bin.points)
        {
            const auto element{ prepared.elements[point] };
            const auto spot_index{ prepared.configured_spot_indices[point] };
            const auto label{ spot_index ? std::string(kConfiguredUmapSpots[*spot_index].label)
                : GetOtherUmapSpotLabel(element) };
            ++spot_counts[label];
            if (include_element_ratios)
            {
                carbon_count += element == Element::CARBON;
                nitrogen_count += element == Element::NITROGEN;
                oxygen_count += element == Element::OXYGEN;
            }
        }
        output << bin.number << ',' << index[0] << ',' << index[1] << ',' << bin.points.size();
        for (const auto feature : { 4u, 5u, 6u })
        {
            const auto & median{ bin.feature_medians[feature] };
            output << ',';
            if (median) output << *median;
        }
        std::vector<std::pair<std::string_view, std::size_t>> ranked_spots(
            spot_counts.begin(), spot_counts.end());
        std::sort(ranked_spots.begin(), ranked_spots.end(), [](const auto & lhs, const auto & rhs)
        {
            if (lhs.second != rhs.second) return lhs.second > rhs.second;
            return lhs.first < rhs.first;
        });
        output << ',' << ranked_spots.size();
        for (std::size_t rank = 0; rank < 3; ++rank)
        {
            if (rank >= ranked_spots.size())
            {
                output << ",,0,0";
                continue;
            }
            const auto & [label, count]{ ranked_spots[rank] };
            output << ',' << label << ',' << count << ','
                << static_cast<double>(count) / static_cast<double>(bin.points.size());
        }
        if (include_element_ratios)
        {
            const double point_count{ static_cast<double>(bin.points.size()) };
            output << ',' << static_cast<double>(carbon_count) / point_count
                << ',' << static_cast<double>(nitrogen_count) / point_count
                << ',' << static_cast<double>(oxygen_count) / point_count;
        }
        output << '\n';
    }
    output.close();
    if (!output)
    {
        error_message = "Failed while writing UMAP mean shift bin output file '" + output_path.string() + "'.";
        return false;
    }
    return true;
}

} // namespace

std::filesystem::path GetEmbeddingOutputPath(
    const std::filesystem::path & output_dir,
    const std::string & model_key_tag)
{
    return output_dir / ("umap_embedding_" + path_helper::EnsureSanitizedTag(model_key_tag) + ".csv");
}

bool WriteEmbeddingOutputs(
    const std::filesystem::path & output_path,
    const PreparedUmapInput & prepared,
    const UmapEmbeddingResult & first,
    const UmapEmbeddingResult & second,
    std::string & error_message)
{
    if (!WriteEmbedding(output_path, prepared, first, second, error_message)) return false;
    Logger::Log(LogLevel::Info, "UMAP embedding written to: " + output_path.string());

    const auto output_stem{ output_path.parent_path() / output_path.stem() };
    if (!WriteMeanShiftBins(prepared, first.binning,
        output_stem.string() + "_mean_shift_bins.csv", false, error_message))
    {
        return false;
    }
    Logger::Log(LogLevel::Info,
        "UMAP mean shift bin count: " + std::to_string(first.binning.bins.size()));

    if (!WriteMeanShiftBins(prepared, second.binning,
        output_stem.string() + "_second_umap_mean_shift_bins.csv", true, error_message))
    {
        return false;
    }
    Logger::Log(LogLevel::Info,
        "Second UMAP mean shift bin count: " + std::to_string(second.binning.bins.size()));
    return true;
}

} // namespace rhbm_gem::core::classifier
