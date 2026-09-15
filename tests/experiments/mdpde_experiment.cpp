#include "support/MDPDEExperiment.hpp"
#include "support/ForwardModelExperiment.hpp"
#include "support/SolverFailureCapture.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <boost/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>

int main(int argc, char ** argv)
{
    namespace fs = std::filesystem;
    try
    {
        if (argc == 6 && std::string(argv[1]) == "forward")
        {
            second_stage_test::RunForwardExperiment(argv[2], argv[3], argv[4], argv[5]);
            return 0;
        }
        if (argc != 4 || std::string(argv[1]) != "solve")
            throw std::runtime_error("Usage: mdpde_experiment solve CAPTURES OUTPUT | forward MANIFEST MAP CAPTURES OUTPUT");
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
