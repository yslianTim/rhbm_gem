#include "support/JointPartialSelection.hpp"
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
#include <memory>
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
    std::thread m_thread;

    void Observe()
    {
        const auto rss = n::CurrentRssBytes();
        if (!rss) return;
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

    std::size_t Peak()
    {
        Observe();
        return m_peak.load();
    }
};

joint_partial_test::Fixture Make(const std::string & name)
{
    joint_partial_test::Fixture f;
    const int components = name == "multi" ? 3 : 1, count = name == "halo" ? 5 : 3;
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    for (int group = 0; group < components; ++group)
        for (int k = 0; k < count; ++k)
        {
            auto atom = std::make_unique<rg::AtomObject>();
            atom->SetSerialID(group * count + k + 1);
            atom->SetElement(Element::CARBON);
            atom->SetPosition(12.0 * group + 0.7 * k, 0, 0);
            atom->SetChainID("A");
            atom->SetComponentID("ALA");
            atom->SetAtomID("CA");
            atoms.push_back(std::move(atom));
            f.a.push_back(2 + 0.2 * k);
            f.b.push_back(0.5 + 0.03 * k);
            f.c.push_back(0.1 - 0.02 * k);
        }
    f.model = std::make_unique<rg::ModelObject>(std::move(atoms));
    f.model->SelectAtoms([&](const auto & atom) {
        return name == "full" || (atom.GetSerialID() - 1) % count == 0;
    });
    const std::array<int, 3> dims{
        static_cast<int>((12.0 * (components - 1) + 0.7 * (count - 1) + 6.4) / 0.16) + 1, 41, 41};
    const std::array<double, 3> spacing{0.16, 0.16, 0.16}, origin{-3.2, -3.2, -3.2};
    auto values = std::make_unique<double[]>(static_cast<std::size_t>(dims[0] * dims[1] * dims[2]));
    for (int z = 0; z < dims[2]; ++z)
        for (int y = 0; y < dims[1]; ++y)
            for (int x = 0; x < dims[0]; ++x)
            {
                const auto position = c::simulation::GridPosition({x, y, z}, spacing, origin);
                double value = 0.0001 * std::sin(0.7 * x + 0.3 * y + 0.2 * z);
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
        if (argc != 4) throw std::invalid_argument(
            "Usage: joint_postprocessing_benchmark full|halo|multi workflow|post DATABASE");
        const std::string name = argv[1], mode = argv[2];
        if ((name != "full" && name != "halo" && name != "multi") ||
            (mode != "workflow" && mode != "post")) return 2;

        auto f = Make(name);
        c::FitOptions options;
        options.estimator = c::PotentialEstimator::JOINT_COMPONENTS;
        options.quiet_mode = true;
        options.thread_size = 1;

        const auto start = Clock::now();
        RssMonitor pipeline_rss;
        c::RunPotentialFittingWorkflow(*f.map, *f.model, options);
        const auto pipeline_seconds = Seconds(start);
        const auto pipeline_peak_rss_bytes = pipeline_rss.Peak();

        const auto save_start = Clock::now();
        { rg::DataRepository repository(argv[3]); repository.SaveModel(*f.model, "model"); }
        const auto save_seconds = Seconds(save_start);

        rusage usage{};
        getrusage(RUSAGE_SELF, &usage);
#ifdef __APPLE__
        const auto process_peak_rss_bytes = static_cast<std::size_t>(usage.ru_maxrss);
#else
        const auto process_peak_rss_bytes = static_cast<std::size_t>(usage.ru_maxrss) * 1024;
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
            if (local)
            {
                const auto & model = local->fit.mdpde.GetModel();
                local_point = {model.GetAmplitude(), model.GetWidth()};
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
                {"joint_point", joint_point},
                {"local_point", local_point},
                {"local_information", information},
                {"local_status", local ? static_cast<int>(local->status) : -1},
                {"local_reason", local ? local->reason : "missing-local-result"}});
        }

        j::object phases{
            {"post_joint_pipeline_seconds", pipeline_seconds},
            {"post_joint_pipeline_peak_rss_bytes", pipeline_peak_rss_bytes},
            {"save_seconds", save_seconds},
            {"total_seconds", Seconds(start)}};
        std::cout << j::serialize(j::object{
            {"case", name},
            {"mode", mode},
            {"phases", phases},
            {"process_peak_rss_bytes", process_peak_rss_bytes},
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
