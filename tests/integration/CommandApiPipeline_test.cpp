#include "utils/domain/FileFingerprint.hpp"
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <sstream>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/CommandSystem.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rg = rhbm_gem;
namespace rgc = rhbm_gem::core;

namespace {

std::filesystem::path FindGeneratedMap(const std::filesystem::path & directory)
{
    for (const auto & entry : std::filesystem::directory_iterator(directory))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".map")
        {
            return entry.path();
        }
    }
    return {};
}

std::size_t CountRegularFiles(const std::filesystem::path & directory)
{
    if (!std::filesystem::exists(directory)) return 0;

    std::size_t count{ 0 };
    for (const auto & entry : std::filesystem::directory_iterator(directory))
    {
        if (entry.is_regular_file())
        {
            ++count;
        }
    }
    return count;
}

void ExpectSelectedAtomsHaveFiniteNonNegativeAlphaR(const rg::ModelObject & model)
{
    const auto & atom_list{ model.GetSelectedAtoms() };
    ASSERT_GT(atom_list.size(), 0u);
    for (const auto * atom : atom_list)
    {
        const auto alpha_r{
            rg::AtomLocalPotentialView::For(*atom).GetAlphaR(
                rg::FittingStage::Second)
        };
        EXPECT_TRUE(std::isfinite(alpha_r));
        EXPECT_GE(alpha_r, 0.0);
    }
}

} // namespace

TEST(CommandApiPipelineTest, FittingIgnoresSimulationTruthAndResolutionMetadata)
{
    command_test::ScopedTempDir directory{"production_truth_isolation"};
    rgc::MapSimulationRequest simulation;
    simulation.model_file_path = command_test::TestDataPath("test_model.cif");
    simulation.output_dir = directory.path();
    simulation.potential_model_choice = rgc::PotentialModel::SINGLE_GAUS;
    simulation.blurring_width_list = {0.5};
    simulation.grid_spacing = 0.2;
    simulation.verbosity = 0;
    ASSERT_TRUE(rgc::RunCommand(simulation).succeeded);
    const auto map{FindGeneratedMap(directory.path())};
    const auto manifest{std::filesystem::path(map.string()+".simulation.json")};
    std::vector<double> reference;
    std::string reference_summary;
    for (int variant = 0; variant < 4; ++variant)
    {
        if (variant == 1) std::filesystem::remove(manifest);
        if (variant == 2) { std::ofstream out(manifest); out << "{broken"; }
        if (variant == 3) { std::ofstream out(manifest); out << R"({"charge_used":99999,"blurring_width":123.456})"; }
        rgc::PotentialAnalysisRequest request;
        request.database_path = directory.path() / (std::to_string(variant)+".sqlite");
        request.model_file_path = simulation.model_file_path;
        request.map_file_path = map;
        request.simulation_flag = true;
        request.simulated_map_resolution = variant == 3 ? 123.456 : 0.5;
        request.saved_key_tag = "isolation";
        request.verbosity = 3;
        testing::internal::CaptureStdout();
        const auto result{rgc::RunCommand(request)};
        const auto log{testing::internal::GetCapturedStdout()};
        ASSERT_TRUE(result.succeeded);
        rg::DataRepository repository{request.database_path};
        auto model{repository.LoadModel("isolation")};
        ASSERT_NE(model, nullptr);
        std::vector<double> values;
        for (const auto & atom : model->GetAtomList())
        {
            const auto view{rg::AtomLocalPotentialView::For(*atom)};
            const auto fit{view.GetGaussianResult(rg::FittingStage::Second)};
            for (const auto & gaussian : {fit.ols.GetModel(), fit.mdpde.GetModel()})
            {
                values.push_back(gaussian.GetAmplitude()); values.push_back(gaussian.GetWidth());
                values.push_back(gaussian.GetOffset());
            }
            values.push_back(view.GetAlphaR(rg::FittingStage::Second));
            for (const auto & samples : {view.GetRawSamplingEntries(false), view.GetPeelingSamplingEntries(false)})
                for (const auto & sample : samples) { values.push_back(sample.point.distance); values.push_back(sample.response); }
        }
        std::istringstream lines{log}; std::string line, summary;
        while (std::getline(lines, line))
            for (const auto * field : {"accepted_iterations", "best_iteration", "stop_reason", "final_uses_polish", "final_state_source"})
                if (line.find(std::string("- ")+field+" =") != std::string::npos) summary += line + '\n';
        ASSERT_FALSE(summary.empty());
        if (variant == 0) { reference = values; reference_summary = summary; }
        else { EXPECT_EQ(values, reference); EXPECT_EQ(summary, reference_summary); }
    }
}

TEST(CommandApiPipelineTest, ExecutesSimulationAnalysisAndDumpPipeline)
{
    command_test::ScopedTempDir temp_dir{ "command_executor_pipeline" };
    const auto maps_dir{ temp_dir.path() / "maps" };
    const auto analysis_output_dir{ temp_dir.path() / "analysis_output" };
    const auto dump_output_dir{ temp_dir.path() / "dump_output" };
    const auto database_path{ temp_dir.path() / "pipeline.sqlite" };

    std::filesystem::create_directories(maps_dir);
    std::filesystem::create_directories(analysis_output_dir);
    std::filesystem::create_directories(dump_output_dir);

    rgc::MapSimulationRequest simulation_request;
    simulation_request.output_dir = maps_dir;
    simulation_request.model_file_path = command_test::TestDataPath("test_model.cif");
    simulation_request.map_file_name = "sim_map";
    simulation_request.blurring_width_list = { 1.50 };

    const auto simulation_result{
        rgc::RunCommand(simulation_request)
    };
    ASSERT_TRUE(simulation_result.succeeded);

    const auto generated_map_file{ FindGeneratedMap(maps_dir) };
    ASSERT_FALSE(generated_map_file.empty());
    ASSERT_TRUE(std::filesystem::exists(generated_map_file));

    rgc::PotentialAnalysisRequest analysis_request;
    analysis_request.database_path = database_path;
    analysis_request.output_dir = analysis_output_dir;
    analysis_request.model_file_path = command_test::TestDataPath("test_model.cif");
    analysis_request.map_file_path = generated_map_file;
    analysis_request.saved_key_tag = "pipeline/model test";

    const auto analysis_result{
        rgc::RunCommand(analysis_request)
    };
    ASSERT_TRUE(analysis_result.succeeded);
    EXPECT_FALSE(std::filesystem::exists(
        analysis_output_dir / "local_fitting_result_pipeline_model_test.csv"));

    rg::DataRepository repository{ database_path };
    auto model{ repository.LoadModel("pipeline/model test") };
    ASSERT_NE(model, nullptr);
    ExpectSelectedAtomsHaveFiniteNonNegativeAlphaR(*model);
    ASSERT_EQ(model->GetAtomList().size(), 1u);
    EXPECT_TRUE(std::isfinite(model->GetAtomList().front()->GetStandardQScore()));
    EXPECT_DOUBLE_EQ(
        model->GetStandardAverageQScore(),
        model->GetAtomList().front()->GetStandardQScore());

    rgc::ResultDumpRequest dump_request;
    dump_request.database_path = database_path;
    dump_request.output_dir = dump_output_dir;
    dump_request.printer_choice = rgc::PrinterType::GAUS_ESTIMATES;
    dump_request.model_key_tag_list = { "pipeline/model test" };

    const auto dump_result{
        rgc::RunCommand(dump_request)
    };
    ASSERT_TRUE(dump_result.succeeded);
    EXPECT_GT(CountRegularFiles(dump_output_dir), 0u);

    const auto analysis_view{ model->GetAnalysisView() };
    const auto group_keys{ analysis_view.CollectAtomGroupKeys() };
    ASSERT_EQ(group_keys.size(), 1u);
    const auto group_key{ group_keys.front() };
    rg::GroupGaussianResult group_result;
    group_result.alpha_g = analysis_view.GetAtomAlphaG(group_key);
    group_result.mean = analysis_view.GetAtomGroupMean(group_key);
    group_result.mdpde = analysis_view.GetAtomGroupMDPDE(group_key);
    group_result.prior = analysis_view.GetAtomGroupPriorWithUncertainty(group_key);
    for (const auto * atom : analysis_view.GetAtomObjectList(group_key))
    {
        const auto & member{ rg::AtomLocalPotentialView::For(*atom).GetGroupMemberResult() };
        ASSERT_TRUE(member.has_value());
        group_result.member_results.emplace_back(*member);
    }
    group_result.member_results.front().is_outlier = true;
    model->EditAnalysis().ApplyAtomGroupGaussianResult(group_key, group_result);
    repository.SaveModel(*model, "pipeline/model test");
    dump_request.printer_choice = rgc::PrinterType::ATOM_OUTLIER;
    ASSERT_TRUE(rgc::RunCommand(dump_request).succeeded);

    bool found_outlier_csv{ false };
    for (const auto & entry : std::filesystem::directory_iterator(dump_output_dir))
    {
        if (entry.path().extension() != ".csv"
            || !entry.path().filename().string().starts_with("atom_outlier_list_")) continue;
        found_outlier_csv = true;
        std::ifstream file{ entry.path() };
        std::string line;
        ASSERT_TRUE(static_cast<bool>(std::getline(file, line)));
        EXPECT_EQ(line, "SerialID,Residue,Element,Spot");
        ASSERT_TRUE(static_cast<bool>(std::getline(file, line)));
        EXPECT_TRUE(line.starts_with(std::to_string(
            analysis_view.GetAtomObjectList(group_key).front()->GetSerialID()) + ","));
    }
    EXPECT_TRUE(found_outlier_csv);

}

TEST(CommandApiPipelineTest, JointOptInSavesTheDirectEndpointAndExportsWithoutSourceInputs)
{
    command_test::ScopedTempDir directory{"joint_command_pipeline"};
    rgc::MapSimulationRequest simulation;
    simulation.model_file_path=directory.path()/"input.cif";
    std::filesystem::copy_file(command_test::TestDataPath("test_model.cif"),simulation.model_file_path);
    {
        std::ifstream original(simulation.model_file_path);
        std::string text((std::istreambuf_iterator<char>(original)),{});
        text.insert(text.rfind('#'),"ATOM 2 C CB . ALA A 1 1.2 0.0 0.0 1.0 0.0 1\n");
        std::ofstream model(simulation.model_file_path); model << text;
    }
    simulation.output_dir=directory.path(); simulation.grid_spacing=.3;
    simulation.potential_model_choice=rgc::PotentialModel::SINGLE_GAUS;
    simulation.blurring_width_list={.5}; simulation.verbosity=0;
    ASSERT_TRUE(rgc::RunCommand(simulation).succeeded);
    const auto map_path=FindGeneratedMap(directory.path()); ASSERT_FALSE(map_path.empty());
    rgc::PotentialAnalysisRequest request;
    request.model_file_path=simulation.model_file_path; request.map_file_path=map_path;
    request.estimator=rgc::PotentialEstimator::JOINT_COMPONENTS; request.only_backbone=true;
    request.database_path=directory.path()/"joint.sqlite"; request.saved_key_tag="joint/model";
    request.verbosity=0;
    for(bool normalization:{false,true})
    {
        request.map_normalization_flag=normalization;
        ASSERT_TRUE(rgc::RunCommand(request).succeeded);
        rg::DataRepository repository{request.database_path};
        auto loaded=repository.LoadModel(request.saved_key_tag);
        ASSERT_TRUE(loaded->GetAnalysisView().GetJointResult());
        const auto & saved=*loaded->GetAnalysisView().GetJointResult();
        ASSERT_TRUE(saved.metadata.map_normalization); ASSERT_TRUE(saved.metadata.software);
        ASSERT_TRUE(saved.selection_domain); EXPECT_EQ(saved.selection_domain->target_indices,(std::vector<std::size_t>{0}));
        EXPECT_EQ(saved.atom_ids,(std::vector<std::string>{"1","2"}));
        EXPECT_EQ(loaded->GetSelectedAtomCount(),1);
        EXPECT_EQ(saved.metadata.map_normalization->requested,normalization);
        EXPECT_EQ(saved.metadata.map_normalization->applied,normalization);
        EXPECT_EQ(saved.metadata.model_sha256,rg::FileSha256(request.model_file_path));
        EXPECT_EQ(saved.metadata.map_sha256,rg::FileSha256(map_path));
        const double sd=rg::ReadMap(map_path)->GetMapValueSD();
        EXPECT_DOUBLE_EQ(saved.metadata.map_normalization->divisor,normalization ? sd : 1);
        auto map=rg::ReadMap(map_path); auto model=rg::ReadModel(request.model_file_path);
        model->SelectAllAtoms(); model->ApplyBackboneSelection(true); if(normalization) map->MapValueArrayNormalization();
        const auto direct=rgc::EstimateJointComponents(*map,*model);
        auto expected=rgc::CaptureJointAnalysisResult(direct,saved.metadata); expected.costs=saved.costs;
        EXPECT_EQ(rg::joint_result_io::Encode(expected),rg::joint_result_io::Encode(saved));
    }
    for(int option=0;option<3;++option)
    {
        auto invalid=request;
        if(option==0) invalid.only_backbone=true;
        if(option==1) invalid.asymmetry_flag=true;
        if(option==2) invalid.sampling_method=SphereSamplingMethod::VolumeUniformRandom;
        EXPECT_EQ(rgc::RunCommand(invalid).succeeded,option!=2);
    }
    {
        auto zero_map=rg::ReadMap(map_path);
        zero_map->SetMapValueArray(std::make_unique<double[]>(zero_map->GetMapValueArraySize()));
        const auto zero_path=directory.path()/"zero.mrc"; rg::WriteMap(zero_path,*zero_map);
        auto zero_request=request; zero_request.map_file_path=zero_path;
        zero_request.map_normalization_flag=true; zero_request.saved_key_tag="zero";
        ASSERT_TRUE(rgc::RunCommand(zero_request).succeeded);
        rg::DataRepository repository{request.database_path};
        auto zero_model=repository.LoadModel("zero");
        ASSERT_TRUE(zero_model->GetAnalysisView().GetJointResult());
        const auto & zero=*zero_model->GetAnalysisView().GetJointResult();
        EXPECT_NE(zero.runtime_convergence,rg::JointCheckStatus::Passed);
        ASSERT_TRUE(zero.metadata.map_normalization);
        EXPECT_TRUE(zero.metadata.map_normalization->requested);
        EXPECT_FALSE(zero.metadata.map_normalization->applied);
        EXPECT_DOUBLE_EQ(zero.metadata.map_normalization->divisor,1);
    }
    {
        auto simulated=request; simulated.simulation_flag=true; simulated.simulated_map_resolution=.5;
        simulated.saved_key_tag="simulation";
        ASSERT_TRUE(rgc::RunCommand(simulated).succeeded);
        rg::DataRepository repository{request.database_path};
        const auto saved=repository.LoadModel("simulation");
        const auto & metadata=saved->GetAnalysisView().GetJointResult()->metadata;
        ASSERT_TRUE(metadata.map_normalization);
        EXPECT_TRUE(metadata.map_normalization->requested);
        EXPECT_FALSE(metadata.map_normalization->applied);
        EXPECT_DOUBLE_EQ(metadata.map_normalization->divisor,1);
    }
    std::filesystem::remove(map_path);
    std::filesystem::remove(request.model_file_path);
    rgc::ResultDumpRequest dump; dump.database_path=request.database_path;
    dump.model_key_tag_list={request.saved_key_tag}; dump.output_dir=directory.path()/"export";
    dump.printer_choice=rgc::PrinterType::JOINT_ESTIMATES; dump.verbosity=0;
    ASSERT_TRUE(rgc::RunCommand(dump).succeeded);
    EXPECT_TRUE(std::filesystem::exists(dump.output_dir/"joint_result_joint_model.json"));
    EXPECT_TRUE(std::filesystem::exists(dump.output_dir/"joint_atoms_joint_model.csv"));
    dump.printer_choice=rgc::PrinterType::GAUS_ESTIMATES; EXPECT_TRUE(rgc::RunCommand(dump).succeeded);
    dump.printer_choice=rgc::PrinterType::ATOM_OUTLIER; EXPECT_TRUE(rgc::RunCommand(dump).succeeded);
    dump.printer_choice=rgc::PrinterType::ATOM_POSITION; EXPECT_TRUE(rgc::RunCommand(dump).succeeded);
    rgc::PotentialDisplayRequest display; display.database_path=request.database_path;
    display.model_key_tag_list={request.saved_key_tag}; display.verbosity=0;
    display.painter_choice=rgc::PainterType::ATOM;
    display.output_dir=directory.path()/"display";
    EXPECT_TRUE(rgc::RunCommand(display).succeeded);
#ifdef HAVE_ROOT
    EXPECT_GT(command_test::CountFilesWithExtension(display.output_dir,".pdf"),0u);
    for (const auto painter : {rgc::PainterType::GAUS,rgc::PainterType::COMPARISON,rgc::PainterType::DEMO})
    {
        display.painter_choice=painter;
        display.verbosity=0;
        display.reference_model_groups={{"with_charge",{request.saved_key_tag}},{"no_charge",{request.saved_key_tag}}};
        EXPECT_TRUE(rgc::RunCommand(display).succeeded) << static_cast<int>(painter);
    }
#endif
#ifdef RHBM_GEM_ENABLE_UMAP
    rgc::UmapEmbeddingRequest umap; umap.database_path=request.database_path;
    umap.model_key_tag=request.saved_key_tag; umap.verbosity=0;
    EXPECT_FALSE(rgc::RunCommand(umap).succeeded);
#endif
    rg::DataRepository repository{request.database_path};
    auto loaded=repository.LoadModel(request.saved_key_tag);
    {
        rgc::JointProblemInput missing_input; missing_input.atom_ids={"1","2"}; missing_input.support.resize(2);
        missing_input.selection_domain=rg::JointSelectionDomain{}; missing_input.selection_domain->target_indices={0};
        auto missing_snapshot=rgc::CaptureJointAnalysisResult(rgc::FitJointComponents(rgc::JointProblem(missing_input),{0.5,0.5}));
        rg::ModelObject missing_model(*loaded); missing_model.EditAnalysis().Clear(); missing_model.EditAnalysis().SetJointResult(missing_snapshot);
        repository.SaveModel(missing_model,"missing");
        display.model_key_tag_list={"missing"}; display.painter_choice=rgc::PainterType::ATOM; display.reference_model_groups.clear();
        EXPECT_TRUE(rgc::RunCommand(display).succeeded);
        dump.model_key_tag_list={"missing"}; dump.printer_choice=rgc::PrinterType::GAUS_ESTIMATES;
        EXPECT_TRUE(rgc::RunCommand(dump).succeeded);
#ifdef RHBM_GEM_ENABLE_UMAP
        umap.model_key_tag="missing"; EXPECT_FALSE(rgc::RunCommand(umap).succeeded);
#endif
    }
    repository.SaveModel(*loaded,"joint_model");
    dump.printer_choice=rgc::PrinterType::JOINT_ESTIMATES;
    dump.model_key_tag_list={request.saved_key_tag,"joint_model"};
    EXPECT_FALSE(rgc::RunCommand(dump).succeeded);
    loaded->EditAnalysis().ClearJointResult(); repository.SaveModel(*loaded,"legacy");
    dump.model_key_tag_list={"legacy"}; EXPECT_FALSE(rgc::RunCommand(dump).succeeded);
}
