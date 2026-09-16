#include "support/MDPDEExperiment.hpp"
#include "support/ForwardModelExperiment.hpp"
#include "support/ObservationMatchedExperiment.hpp"
#include "support/EstimatedNeighborSweep.hpp"
#include "support/SolverFailureCapture.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <boost/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>

int main(int argc, char ** argv)
{
    namespace fs = std::filesystem;
    try
    {
        if (argc == 6 && std::string(argv[1]) == "matched-sweep")
        {
            second_stage_test::matched::RunEstimatedNeighborSweep(argv[2], argv[3], argv[4], argv[5]);
            return 0;
        }
        if (argc == 6 && std::string(argv[1]) == "matched")
        {
            second_stage_test::matched::Run(argv[2], argv[3], argv[4], argv[5]);
            return 0;
        }
        if (argc == 6 && std::string(argv[1]) == "forward")
        {
            second_stage_test::RunForwardExperiment(argv[2], argv[3], argv[4], argv[5]);
            return 0;
        }
        if (argc == 5 && std::string(argv[1]) == "refine")
        {
            const int budget{std::stoi(argv[4])};
            const fs::path output{argv[3]}; fs::create_directories(output);
            boost::json::array index;
            auto run = [&](second_stage_test::ShapeFixture fixture, const std::string & name, bool replay) {
                fixture.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(fixture.alpha,fixture.dataset,fixture.options);
                auto result{second_stage_test::RefineMDPDEEndpoint(fixture,budget)};
                result.evidence["fixture"] = name; result.evidence["exact_replay"] = replay;
                const auto target{output / (name + ".json")};
                std::ofstream out(target); out.exceptions(std::ios::failbit | std::ios::badbit);
                out << boost::json::serialize(result.evidence); index.emplace_back(target.filename().string());
            };
            std::size_t captures{};
            for (const auto & file : fs::directory_iterator(argv[2]))
            {
                if (file.path().extension() != ".txt" || !file.path().filename().string().starts_with("shape-")) continue;
                if (!second_stage_test::ReplaySolverFailure(file.path().string()))
                    throw std::runtime_error("Exact replay failed: " + file.path().string());
                run(second_stage_test::ReadShapeFixture(file.path().string()),file.path().stem().string(),true); ++captures;
            }
            if (!captures) throw std::runtime_error("No refinement captures.");
            for (int k = 0; k < 4; ++k)
            {
                second_stage_test::ShapeFixture fixture;
                fixture.alpha = std::array{0.0,0.1,0.2,0.5}[static_cast<std::size_t>(k)];
                fixture.dataset.X.resize(40,2); fixture.dataset.y.resize(40);
                for (int i = 0; i < 40; ++i)
                {
                    const double r{static_cast<double>(i)/40.0};
                    fixture.dataset.X.row(i) << 1.0,-0.5*r*r;
                    fixture.dataset.y(i) = 1.2-2.0*r*r+(0.003+0.005*k)*std::sin(3.0*i+0.2+k);
                }
                run(std::move(fixture),"independent-"+std::to_string(k),false);
            }
            std::ofstream(output / "index.json") << boost::json::serialize(boost::json::object{
                {"schema_version",1},{"shape_count",captures},{"independent_count",4},{"files",index}});
            return 0;
        }
        if (argc != 4 || std::string(argv[1]) != "solve")
            throw std::runtime_error("Usage: mdpde_experiment solve CAPTURES OUTPUT | forward|matched MANIFEST MAP CAPTURES OUTPUT | matched-sweep MANIFEST MAP STATE_INDEX OUTPUT | refine CAPTURES OUTPUT BUDGET");
        const fs::path output{ argv[3] };
        fs::create_directories(output);
        std::size_t count{}, offsets{};
        boost::json::array index;
        for (const auto & file : fs::directory_iterator(argv[2]))
        {
            if (file.path().extension() != ".txt") continue;
            if (!second_stage_test::ReplaySolverFailure(file.path().string()))
                throw std::runtime_error("Exact replay failed: " + file.path().string());
            if (file.path().filename().string().starts_with("offset-")) { ++offsets; continue; }
            auto result{ second_stage_test::CompareMDPDE(second_stage_test::ReadShapeFixture(file.path().string())) };
            result["fixture"] = file.path().filename().string();
            const auto target{ output / (file.path().stem().string() + ".json") };
            std::ofstream out(target); out.exceptions(std::ios::failbit | std::ios::badbit);
            out << boost::json::serialize(result);
            index.emplace_back(target.filename().string()); ++count;
        }
        if (!count) throw std::runtime_error("No shape fixtures; external evidence cannot be skipped.");
        std::ofstream(output / "index.json") << boost::json::serialize(boost::json::object{
            {"schema_version", 1}, {"shape_count", count}, {"offset_exact_replays", offsets}, {"files", index}});
        std::cout << "Compared " << count << " shape fixtures; replayed " << offsets << " offset fixtures.\n";
        return 0;
    }
    catch (const std::exception & e) { std::cerr << e.what() << '\n'; return 1; }
}
