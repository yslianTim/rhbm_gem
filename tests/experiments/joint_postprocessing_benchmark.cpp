#include "support/JointPartialSelection.hpp"
#include "core/detail/PotentialFittingWorkflow.hpp"
#include "core/detail/joint_component/ResourceWork.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <boost/json.hpp>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <thread>
#include <vector>

namespace {
namespace rg = rhbm_gem;
namespace c = rg::core;
namespace n = c::joint_component;
namespace j = boost::json;
using Clock = std::chrono::steady_clock;

double Seconds(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

class RssMonitor
{
    std::atomic<bool> m_running{true};
    std::atomic<std::size_t> m_peak{};
    std::atomic<bool> m_available{};
    std::thread m_thread;

    void Observe()
    {
        const auto rss = n::CurrentRssBytes();
        if (!rss) return;
        m_available = true;
        auto peak = m_peak.load();
        while (peak < *rss && !m_peak.compare_exchange_weak(peak, *rss)) {}
    }

public:
    RssMonitor() : m_thread([this] {
        while (m_running.load())
        {
            Observe();
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        Observe();
    }) { Observe(); }

    ~RssMonitor()
    {
        m_running = false;
        m_thread.join();
    }

    std::optional<std::size_t> Peak()
    {
        Observe();
        return m_available.load() ? std::optional<std::size_t>{m_peak.load()} : std::nullopt;
    }
};

joint_partial_test::Fixture Make(const std::string & name)
{
    joint_partial_test::Fixture f;
    const int components = name == "multi" ? 3 : 1;
    const int count = name == "halo" ? 5 : name == "target-heavy" ? 6 : name == "connected" ? 9 : 3;
    const bool full_selection = name != "halo";
    const double grid_spacing = name == "dense" ? 0.08 : 0.16;
    const double atom_spacing = name == "target-heavy" ? 0.28 : 0.7;
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    for (int group = 0; group < components; ++group)
        for (int k = 0; k < count; ++k)
        {
            auto atom = std::make_unique<rg::AtomObject>();
            atom->SetSerialID(group * count + k + 1);
            atom->SetElement(Element::CARBON);
            const double x = name == "target-heavy" ? 0.7 * static_cast<double>(k % 3) : atom_spacing * k;
            const double y = name == "target-heavy" ? 0.7 * static_cast<double>(k / 3) : 0.0;
            atom->SetPosition(12.0 * group + x, y, 0);
            atom->SetChainID("A");
            atom->SetComponentID("ALA");
            atom->SetAtomID("CA");
            atoms.push_back(std::move(atom));
            f.a.push_back(2 + 0.2 * k);
            f.b.push_back(0.5 + 0.03 * k);
            f.c.push_back(0.1 - 0.02 * k);
        }
    f.model = std::make_unique<rg::ModelObject>(std::move(atoms));
    f.model->SelectAtoms([&](const auto & atom) { return full_selection || (atom.GetSerialID() - 1) % count == 0; });
    const std::array<int, 3> dims{
        static_cast<int>(std::ceil((12.0 * (components - 1) + atom_spacing * (count - 1) + 6.4) / grid_spacing)) + 1,
        static_cast<int>(std::ceil(6.4 / grid_spacing)) + 1,
        static_cast<int>(std::ceil(6.4 / grid_spacing)) + 1};
    const std::array<double, 3> spacing{grid_spacing, grid_spacing, grid_spacing}, origin{-3.2, -3.2, -3.2};
    auto values = std::make_unique<double[]>(static_cast<std::size_t>(dims[0] * dims[1] * dims[2]));
    for (int z = 0; z < dims[2]; ++z)
        for (int y = 0; y < dims[1]; ++y)
            for (int x = 0; x < dims[0]; ++x)
            {
                const auto position = c::simulation::GridPosition({x, y, z}, spacing, origin);
                const double noise_scale = name == "noisy" ? 0.002 : 0.0001;
                double value = noise_scale * std::sin(0.7 * x + 0.3 * y + 0.2 * z);
                if (name == "outlier" && x == 20 && y == 20 && z == 20) value += 0.1;
                for (std::size_t atom = 0; atom < f.a.size(); ++atom)
                {
                    const auto kernel = c::joint_component::EvaluateKernel(
                        c::simulation::SupportSquare(position,
                            f.model->GetAtomList()[atom]->GetPosition()), f.b[atom], 2.5);
                    value += f.a[atom] * kernel.gaussian + f.c[atom] * kernel.charge;
                }
                values[static_cast<std::size_t>(x + dims[0] * (y + dims[1] * z))] = value;
            }
    f.map = std::make_unique<rg::MapObject>(dims, spacing, origin, std::move(values));
    return f;
}
}

int main(int argc, char ** argv)
{
    try
    {
        if (argc != 3) throw std::invalid_argument(
            "Usage: joint_postprocessing_benchmark full|halo|multi|noisy|outlier|dense|target-heavy|connected DATABASE");
        const std::string name = argv[1];
        if (name != "full" && name != "halo" && name != "multi" && name != "noisy" &&
            name != "outlier" && name != "dense" && name != "target-heavy" && name != "connected") return 2;
        const std::string policy = "pure";

        auto f = Make(name);
        c::FitOptions options;
        options.estimator = c::PotentialEstimator::JOINT_COMPONENTS;
        options.quiet_mode = true;
        options.thread_size = 1;

        const auto total_start = Clock::now();
        c::detail::JointPostprocessingMetrics postprocessing_metrics;
        double pipeline_seconds{};
        std::optional<std::size_t> pipeline_peak_rss_bytes;
        {
            const auto workflow_start = Clock::now();
            RssMonitor pipeline_rss;
            c::detail::RunJointComponentWorkflow(*f.map, *f.model, options, &postprocessing_metrics);
            pipeline_seconds = Seconds(workflow_start);
            pipeline_peak_rss_bytes = pipeline_rss.Peak();
        }

        const auto save_start = Clock::now();
        { rg::DataRepository repository(argv[2]); repository.SaveModel(*f.model, "model"); }
        const auto save_seconds = Seconds(save_start);

        rusage usage{};
        const bool process_peak_available = getrusage(RUSAGE_SELF, &usage) == 0;
#ifdef __APPLE__
        const auto process_peak_rss_bytes = process_peak_available ?
            std::optional<std::size_t>{static_cast<std::size_t>(usage.ru_maxrss)} : std::nullopt;
#else
        const auto process_peak_rss_bytes = process_peak_available ?
            std::optional<std::size_t>{static_cast<std::size_t>(usage.ru_maxrss) * 1024} : std::nullopt;
#endif

        j::array endpoints, targets;
        const auto & snapshot = *f.model->GetAnalysisView().GetJointResult();
        for (const auto & component : snapshot.components)
            endpoints.push_back(j::object{
                {"id", component.id},
                {"convergence", static_cast<int>(component.runtime_convergence)},
                {"ac", component.state ? j::value_from(component.state->ac) : j::value{}},
                {"b", component.state ? j::value_from(component.state->b) : j::value{}}});
        for (const auto * atom : f.model->GetSelectedAtoms())
        {
            const auto view = rg::AtomLocalPotentialView::For(*atom);
            const auto & joint = view.GetStageEstimate(rg::FittingStage::Second);
            j::array joint_point, local_point, information;
            if (joint.point)
                joint_point = {joint.point->GetAmplitude(), joint.point->GetWidth(), joint.point->GetOffset()};
            const auto & local = view.GetPostJointLocalResult();
            j::object failure_reasons;
            if (local)
            {
                const auto & model = local->fit.mdpde.GetModel();
                local_point = {model.GetAmplitude(), model.GetWidth()};
                std::map<std::string, std::size_t> reason_counts;
                for (const auto & failure : local->sample_failures) ++reason_counts[failure.reason];
                for (const auto & [reason, count] : reason_counts) failure_reasons[reason] = count;
                if (local->information)
                {
                    for (int r = 0; r < 2; ++r)
                        for (int col = 0; col < 2; ++col)
                            information.push_back(local->information->precision(r, col));
                    for (int row = 0; row < 2; ++row)
                        information.push_back(local->information->moment(row));
                }
            }
            targets.push_back(j::object{
                {"id", atom->GetSerialID()},
                {"joint_convergence", static_cast<int>(joint.convergence)},
                {"joint_point", joint_point},
                {"local_point", local_point},
                {"local_information", information},
                {"local_status", local ? static_cast<int>(local->status) : -1},
                {"local_reason", local ? local->reason : "missing-local-result"},
                {"local_sample_failure_reasons", failure_reasons}});
        }

        j::object phases{
            {"postprocessing_seconds", postprocessing_metrics.total_postprocessing_seconds},
            {"context_construction_seconds", postprocessing_metrics.context_construction_seconds},
            {"peeling_seconds", postprocessing_metrics.peeling_seconds},
            {"charge_correction_seconds", postprocessing_metrics.charge_correction_seconds},
            {"local_alpha_training_seconds", postprocessing_metrics.local_alpha_training_seconds},
            {"local_mdpde_fitting_seconds", postprocessing_metrics.local_mdpde_fitting_seconds},
            {"group_alpha_training_seconds", postprocessing_metrics.group_alpha_training_seconds},
            {"group_web_seconds", postprocessing_metrics.group_web_seconds},
            {"workflow_seconds", pipeline_seconds},
            {"workflow_peak_rss_bytes", pipeline_peak_rss_bytes ? j::value(*pipeline_peak_rss_bytes) : j::value{}},
            {"workflow_peak_rss_available", pipeline_peak_rss_bytes.has_value()},
            {"save_seconds", save_seconds},
            {"total_seconds", Seconds(total_start)}};
        const j::object preparation{
            {"observation_count", snapshot.row_ids.size()},
            {"contributor_count", snapshot.atom_ids.size()},
            {"target_count", f.model->GetSelectedAtoms().size()},
            {"snapshot_identity_validation_count", postprocessing_metrics.snapshot_identity_validation_count},
            {"full_row_identity_scan_count", postprocessing_metrics.full_row_identity_scan_count},
            {"row_mapping_build_count", postprocessing_metrics.observation_row_mapping_builds},
            {"contributor_mapping_build_count", postprocessing_metrics.contributor_mapping_builds},
            {"joint_stage_mapping_build_count", postprocessing_metrics.joint_stage_mapping_builds},
            {"target_preparation_count", postprocessing_metrics.target_preparation_count},
            {"prepared_sample_count", postprocessing_metrics.total_prepared_sample_count},
            {"valid_local_sample_count", postprocessing_metrics.valid_local_sample_count},
            {"invalid_local_sample_count", postprocessing_metrics.invalid_local_sample_count},
            {"raw_peeling_sample_count", postprocessing_metrics.raw_peeling_sample_count},
            {"nonzero_stencil_node_count", postprocessing_metrics.nonzero_stencil_node_count},
            {"available_peeling_sample_count", postprocessing_metrics.available_peeling_sample_count},
            {"unavailable_peeling_sample_count", postprocessing_metrics.unavailable_peeling_sample_count},
            {"map_geometry_boundary_samples", postprocessing_metrics.map_geometry_boundary_samples},
            {"map_geometry_boundary_nodes", postprocessing_metrics.map_geometry_boundary_nodes},
            {"original_support_boundary_samples", postprocessing_metrics.original_support_boundary_samples},
            {"original_support_boundary_nodes", postprocessing_metrics.original_support_boundary_nodes},
            {"fitting_mask_only_intersection_samples", postprocessing_metrics.fitting_mask_only_intersection_samples},
            {"fitting_mask_only_intersection_nodes", postprocessing_metrics.fitting_mask_only_intersection_nodes},
            {"mapping_mismatch_samples", postprocessing_metrics.mapping_mismatch_samples},
            {"missing_contributor_samples", postprocessing_metrics.missing_contributor_samples},
            {"other_peeling_failure_samples", postprocessing_metrics.other_peeling_failure_samples}};
        const j::object memory{
            {"estimated_context_lookup_bytes", postprocessing_metrics.estimated_context_lookup_bytes},
            {"estimated_peak_workspace_payload_bytes", postprocessing_metrics.estimated_peak_postprocessing_payload_bytes},
            {"current_rss_before_postprocessing_bytes", postprocessing_metrics.rss_before_postprocessing_bytes ?
                j::value(*postprocessing_metrics.rss_before_postprocessing_bytes) : j::value{}},
            {"current_rss_after_context_bytes", postprocessing_metrics.rss_after_context_bytes ?
                j::value(*postprocessing_metrics.rss_after_context_bytes) : j::value{}},
            {"current_rss_after_peeling_bytes", postprocessing_metrics.rss_after_peeling_bytes ?
                j::value(*postprocessing_metrics.rss_after_peeling_bytes) : j::value{}},
            {"current_rss_after_local_bytes", postprocessing_metrics.rss_after_local_bytes ?
                j::value(*postprocessing_metrics.rss_after_local_bytes) : j::value{}},
            {"current_rss_after_group_bytes", postprocessing_metrics.rss_after_group_bytes ?
                j::value(*postprocessing_metrics.rss_after_group_bytes) : j::value{}},
            {"process_peak_rss_bytes", process_peak_rss_bytes ? j::value(*process_peak_rss_bytes) : j::value{}},
            {"process_peak_rss_available", process_peak_available},
            {"estimate_scope", "container and sample payload estimate; allocator overhead excluded"},
            {"rss_scope", "process high-water includes fixture and persistence; sampled workflow peak may miss short spikes"}};
        std::cout << j::serialize(j::object{
            {"case", name},
            {"mode", "workflow"},
            {"joint_policy", policy},
            {"phases", phases},
            {"preparation", preparation},
            {"memory", memory},
            {"objective", snapshot.objective ? j::value(*snapshot.objective) : j::value{}},
            {"endpoints", endpoints},
            {"targets", targets}}) << '\n';
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
