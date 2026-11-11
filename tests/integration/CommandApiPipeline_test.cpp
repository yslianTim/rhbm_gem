#include "utils/domain/FileFingerprint.hpp"
#include <gtest/gtest.h>
#include <boost/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/CommandSystem.hpp>
#include "core/estimator/joint_component/JointComponentEstimator.hpp"
#include "core/estimator/PotentialFittingWorkflow.hpp"
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

// Calls the instrumented command entry point in-process. The CLI smoke test runs the production executable.
int RunCli(std::vector<std::string> arguments)
{
    std::vector<char *> argv;
    argv.reserve(arguments.size());
    for (auto & argument : arguments) argv.push_back(argument.data());
    return rgc::RunCommandCLI(static_cast<int>(argv.size()), argv.data());
}

void InjectJointSearchFailure(rg::JointAnalysisResult & fit)
{
    fit.search_completed=false;
    for (auto & component:fit.components)
    {
        component.search_completed=false;
        component.stop_reason="block-recovery-exhausted";
    }
}

void InjectJointNoDonors(rg::JointAnalysisResult & fit)
{
    if (!fit.target_evidence) throw std::runtime_error("Test fixture lacks target evidence.");
    for (auto & point:fit.target_evidence->point_atoms) point.local_conditioning="Warning";
    for (auto & component:fit.components)
        if (component.target_evidence)
            for (auto & point:component.target_evidence->point_atoms)
                point.local_conditioning="Warning";
}

void ExpectJointNumericalResultPreserved(
    const rg::JointAnalysisResult & before,
    const rg::JointAnalysisResult & after)
{
    auto comparable_before=before;
    comparable_before.metadata=after.metadata;
    comparable_before.costs=after.costs;
    comparable_before.halo_median_reconstruction=after.halo_median_reconstruction;
    ASSERT_EQ(comparable_before.components.size(),after.components.size());
    for(std::size_t i=0;i<after.components.size();++i)
        comparable_before.components[i].local_conditioning_seconds=after.components[i].local_conditioning_seconds;
    EXPECT_EQ(rg::joint_result_io::Encode(comparable_before),rg::joint_result_io::Encode(after));
}

struct ScopedJointPostprocessingSnapshotHook
{
    explicit ScopedJointPostprocessingSnapshotHook(rgc::detail::JointPostprocessingSnapshotHook hook)
    {
        rgc::detail::JointPostprocessingSnapshotHookForTesting()=hook;
    }
    ~ScopedJointPostprocessingSnapshotHook()
    {
        rgc::detail::JointPostprocessingSnapshotHookForTesting()=nullptr;
    }
};

void ExpectSelectedAtomsHaveFiniteNonNegativeAlphaR(const rg::ModelObject & model)
{
    const auto & atom_list{ model.GetSelectedAtoms() };
    ASSERT_GT(atom_list.size(), 0u);
    for (const auto * atom : atom_list)
    {
        const auto alpha_r{
            rg::AtomLocalPotentialView::For(*atom).GetInitialStageAlphaR()
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
    std::optional<rg::JointAnalysisResult> reference_joint;
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
        static_cast<void>(testing::internal::GetCapturedStdout());
        ASSERT_TRUE(result.succeeded);
        rg::DataRepository repository{request.database_path};
        auto model{repository.LoadModel("isolation")};
        ASSERT_NE(model, nullptr);
        ASSERT_EQ(result.outcome, rgc::CommandOutcome::JointNotConverged);
        const auto & joint = model->GetAnalysisView().GetJointResult();
        ASSERT_TRUE(joint);
        if (!reference_joint) reference_joint = *joint;
        else ExpectJointNumericalResultPreserved(*reference_joint, *joint);
        std::vector<double> values;
        for (const auto & atom : model->GetAtomList())
        {
            const auto view{rg::AtomLocalPotentialView::For(*atom)};
            const auto fit{view.GetInitialStageGaussianResult()};
            for (const auto & gaussian : {fit.ols.GetModel(), fit.mdpde.GetModel()})
            {
                values.push_back(gaussian.GetAmplitude()); values.push_back(gaussian.GetWidth());
                values.push_back(gaussian.GetOffset());
            }
            values.push_back(view.GetInitialStageAlphaR());
            for (const auto & sample : view.GetRawSamplingEntriesView())
            { values.push_back(sample.point.distance); values.push_back(sample.response); }
        }
        if (variant == 0) reference = values;
        else EXPECT_EQ(values, reference);
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

    ASSERT_TRUE(model->GetAnalysisView().GetJointResult());
    EXPECT_EQ(analysis_result.outcome, rgc::CommandOutcome::JointNotConverged);
    dump_request.printer_choice = rgc::PrinterType::ATOM_OUTLIER;
    EXPECT_TRUE(rgc::RunCommand(dump_request).succeeded);
}

TEST(CommandApiPipelineTest, DirectPotentialFittingExcludesHydrogenAndKeepsSelectedSubset)
{
    command_test::ScopedTempDir directory{ "direct_hydrogen_free_fitting" };
    const auto model_path{ command_test::TestDataPath("potential_analysis_multichain.cif") };
    const auto map_path{ command_test::GenerateMapFile(directory.path() / "maps", model_path) };
    auto map{ rg::ReadMap(map_path) };
    auto model{ rg::ReadModel(model_path) };
    ASSERT_NE(map, nullptr);
    ASSERT_NE(model, nullptr);
    ASSERT_TRUE(map->MapValueArrayNormalization());
    model->SelectAtoms([](const rg::AtomObject & atom)
    {
        return atom.GetSerialID() == 1 || atom.GetSerialID() == 7;
    });

    rgc::FitOptions options;
    options.thread_size = 1;
    options.quiet_mode = true;
    rgc::RunPotentialFittingWorkflow(*map, *model, options);

    ASSERT_EQ(model->GetSelectedAtomCount(), 1u);
    EXPECT_EQ(model->GetSelectedAtoms().front()->GetSerialID(), 1);
    EXPECT_EQ(model->GetSelectedAtoms().front()->GetElement(), Element::CARBON);
    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*model->FindAtomPtr(7)).IsAvailable());
    const auto & joint = model->GetAnalysisView().GetJointResult();
    ASSERT_TRUE(joint);
    ASSERT_TRUE(joint->selection_domain);
    ASSERT_EQ(joint->selection_domain->target_indices.size(), 1u);
    EXPECT_EQ(joint->atom_ids.at(joint->selection_domain->target_indices.front()), "1");
    EXPECT_TRUE(joint->initialization.valid);
    EXPECT_FALSE(rg::IsQualifiedJointAnalysisResult(*joint));

    auto hydrogen_only_model{ rg::ReadModel(model_path) };
    ASSERT_NE(hydrogen_only_model, nullptr);
    hydrogen_only_model->SelectAtoms([](const rg::AtomObject & atom)
    {
        return atom.GetSerialID() == 7;
    });
    EXPECT_THROW(rgc::RunPotentialFittingWorkflow(*map, *hydrogen_only_model, options),
        std::invalid_argument);
    EXPECT_EQ(hydrogen_only_model->GetSelectedAtomCount(), 0u);
}

TEST(CommandApiPipelineTest, JointPersistsSmallSelectionCalibrationFailureAndExportsWithoutSourceInputs)
{
    command_test::ScopedTempDir directory{"joint_command_pipeline"};
    rgc::MapSimulationRequest simulation;
    simulation.model_file_path=directory.path()/"input.cif";
    std::filesystem::copy_file(
        command_test::TestDataPath("potential_analysis_multichain.cif"), simulation.model_file_path);
    simulation.output_dir=directory.path(); simulation.grid_spacing=.3;
    simulation.potential_model_choice=rgc::PotentialModel::SINGLE_GAUS;
    simulation.blurring_width_list={.5}; simulation.verbosity=0;
    ASSERT_TRUE(rgc::RunCommand(simulation).succeeded);
    const auto map_path=FindGeneratedMap(directory.path()); ASSERT_FALSE(map_path.empty());
    rgc::PotentialAnalysisRequest request;
    request.model_file_path=simulation.model_file_path; request.map_file_path=map_path;
    request.selected_chain_ids={"A"};
    request.database_path=directory.path()/"joint.sqlite"; request.saved_key_tag="joint/model";
    request.verbosity=0;
    const auto api_result=rgc::RunCommand(request);
    ASSERT_TRUE(api_result.succeeded);
    rg::DataRepository repository{request.database_path};
    auto loaded=repository.LoadModel(request.saved_key_tag);
    ASSERT_NE(loaded, nullptr);
    ASSERT_EQ(loaded->GetNumberOfAtom(), 8u);
    EXPECT_EQ(loaded->GetSelectedAtomCount(), 2u);
    ASSERT_TRUE(loaded->GetAnalysisView().GetJointResult());
    const auto & saved=*loaded->GetAnalysisView().GetJointResult();
    const auto expected_outcome = !saved.initialization.valid ? rgc::CommandOutcome::JointInitializationFailed :
        !rg::IsQualifiedJointAnalysisResult(saved) ? rgc::CommandOutcome::JointNotConverged :
        rg::HasUnavailableJointHaloMedianReconstruction(saved) ? rgc::CommandOutcome::JointDownstreamUnavailable :
        rgc::CommandOutcome::JointQualified;
    EXPECT_EQ(api_result.outcome,expected_outcome);
    ASSERT_TRUE(saved.metadata.map_normalization); ASSERT_TRUE(saved.metadata.software);
    ASSERT_TRUE(saved.selection_domain);
    EXPECT_EQ(saved.selection_domain->target_indices,(std::vector<std::size_t>{0,1}));
    EXPECT_EQ(saved.atom_ids.front(), "1");
    EXPECT_EQ(saved.atom_ids[1], "2");
    const auto has_contributor = [&saved](const std::string & id)
    {
        return std::find(saved.atom_ids.begin(), saved.atom_ids.end(), id) != saved.atom_ids.end();
    };
    for (const auto * halo_id : {"3", "4", "5", "6"}) EXPECT_TRUE(has_contributor(halo_id));
    EXPECT_FALSE(has_contributor("7"));
    EXPECT_FALSE(has_contributor("8"));
    EXPECT_TRUE(saved.metadata.map_normalization->requested);
    EXPECT_TRUE(saved.metadata.map_normalization->applied);
    EXPECT_EQ(saved.metadata.model_sha256,rg::FileSha256(request.model_file_path));
    EXPECT_EQ(saved.metadata.map_sha256,rg::FileSha256(map_path));
    const double sd=rg::ReadMap(map_path)->GetMapValueSD();
    EXPECT_DOUBLE_EQ(saved.metadata.map_normalization->divisor, sd);
    auto map=rg::ReadMap(map_path); auto model=rg::ReadModel(request.model_file_path);
    ASSERT_NE(map, nullptr); ASSERT_NE(model, nullptr);
    ASSERT_TRUE(map->MapValueArrayNormalization());
    EXPECT_FALSE(rg::IsQualifiedJointAnalysisResult(saved));
    const auto calibration=std::find_if(saved.evidence.begin(),saved.evidence.end(),[](const auto & item) {
        return item.name=="global-width-calibration";
    });
    ASSERT_NE(calibration,saved.evidence.end());
    EXPECT_EQ(calibration->status,rg::JointCheckStatus::Unavailable);
    EXPECT_TRUE(saved.components.empty());
    model->SelectAtoms([](const rg::AtomObject & atom)
    {
        return atom.GetElement() != Element::HYDROGEN && atom.GetChainID() == "A";
    });
    const auto target_problem=rgc::BuildPureHaloSupportExclusionProblem(*map,*model);
    EXPECT_EQ(saved.row_ids,target_problem.OriginalInput().row_ids);
    EXPECT_EQ(target_problem.Input().selection_domain->target_indices,
        (std::vector<std::size_t>{0,1}));

    {
        auto zero_map=rg::ReadMap(map_path);
        zero_map->SetMapValueArray(std::make_unique<double[]>(zero_map->GetMapValueArraySize()));
        const auto zero_path=directory.path()/"zero.mrc"; rg::WriteMap(zero_path,*zero_map);
        auto zero_request=request; zero_request.map_file_path=zero_path;
        zero_request.saved_key_tag="zero";
        const auto zero_result=rgc::RunCommand(zero_request);
        ASSERT_TRUE(zero_result.succeeded);
        EXPECT_EQ(zero_result.outcome,rgc::CommandOutcome::JointInitializationFailed);
        const auto zero_model=repository.LoadModel("zero");
        ASSERT_TRUE(zero_model->GetAnalysisView().GetJointResult());
        EXPECT_FALSE(zero_model->GetAnalysisView().GetJointResult()->initialization.valid);
    }
    {
        auto sparse_map=rg::ReadMap(map_path);
        auto values=std::make_unique<double[]>(sparse_map->GetMapValueArraySize());
        values[sparse_map->GetMapValueArraySize()-1]=1.0;
        sparse_map->SetMapValueArray(std::move(values));
        const auto sparse_path=directory.path()/"sparse.mrc"; rg::WriteMap(sparse_path,*sparse_map);
        auto sparse_request=request; sparse_request.map_file_path=sparse_path;
        sparse_request.saved_key_tag="sparse";
        const auto sparse_result=rgc::RunCommand(sparse_request);
        ASSERT_TRUE(sparse_result.succeeded);
        EXPECT_EQ(sparse_result.outcome,rgc::CommandOutcome::JointInitializationFailed);
        const auto sparse_model=repository.LoadModel("sparse");
        ASSERT_TRUE(sparse_model->GetAnalysisView().GetJointResult());
        EXPECT_FALSE(sparse_model->GetAnalysisView().GetJointResult()->initialization.valid);
        const auto sparse_cli_database=directory.path()/"sparse-cli.sqlite";
        const int sparse_cli_status=RunCli({"RHBM-GEM","potential_analysis",
            "--model",request.model_file_path.string(),"--map",sparse_path.string(),
            "--database",sparse_cli_database.string(),"--save-key","sparse-cli","--select-chain","A",
            "--verbose","0"});
        EXPECT_EQ(sparse_cli_status,2);
    }
    {
        auto simulated=request; simulated.simulation_flag=true; simulated.simulated_map_resolution=.5;
        simulated.saved_key_tag="simulation";
        ASSERT_TRUE(rgc::RunCommand(simulated).succeeded);
        const auto simulation_model=repository.LoadModel("simulation");
        const auto & metadata=simulation_model->GetAnalysisView().GetJointResult()->metadata;
        ASSERT_TRUE(metadata.map_normalization);
        EXPECT_FALSE(metadata.map_normalization->requested);
        EXPECT_FALSE(metadata.map_normalization->applied);
        EXPECT_DOUBLE_EQ(metadata.map_normalization->divisor,1);
    }
    const auto cli_database=directory.path()/"cli.sqlite";
    const auto cli_output=directory.path()/"cli-output";
    const int cli_status=RunCli({"RHBM-GEM","potential_analysis",
        "--model",request.model_file_path.string(),"--map",map_path.string(),
        "--database",cli_database.string(),"--save-key","cli-joint","--select-chain","A",
        "--folder",cli_output.string(),"--verbose","0"});
    rg::DataRepository cli_repository{cli_database};
    const auto cli_model=cli_repository.LoadModel("cli-joint");
    ASSERT_NE(cli_model,nullptr);
    ASSERT_TRUE(cli_model->GetAnalysisView().GetJointResult());
    const auto & cli_fit=*cli_model->GetAnalysisView().GetJointResult();
    const int expected_cli_status=!cli_fit.initialization.valid ? 2 :
        !rg::IsQualifiedJointAnalysisResult(cli_fit) ? 3 :
        rg::HasUnavailableJointHaloMedianReconstruction(cli_fit) ? 4 : 0;
    EXPECT_EQ(cli_status,expected_cli_status);
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
    loaded=repository.LoadModel(request.saved_key_tag);
    {
        const auto missing_snapshot=*loaded->GetAnalysisView().GetJointResult();
        rg::ModelObject missing_model(*loaded); missing_model.EditAnalysis().Clear(); missing_model.EditAnalysis().ApplyJointResult(missing_snapshot,"imported-joint");
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
    loaded->EditAnalysis().ClearJointResult(); repository.SaveModel(*loaded,"no-joint-result");
    dump.model_key_tag_list={"no-joint-result"}; EXPECT_FALSE(rgc::RunCommand(dump).succeeded);
}

TEST(CommandApiPipelineTest, JointFailureAndUnavailableOutcomesPersistThroughApiCliAndReload)
{
    command_test::ScopedTempDir directory{"joint_failure_outcomes"};
    const auto maps_dir=directory.path()/"maps";
    std::filesystem::create_directories(maps_dir);
    const auto model_path=directory.path()/"input.cif";
    const auto template_path=command_test::TestDataPath("test_model.cif");
    std::ifstream template_file(template_path);
    ASSERT_TRUE(template_file);
    const std::string template_text{std::istreambuf_iterator<char>{template_file},
        std::istreambuf_iterator<char>{}};
    const auto atom_marker=template_text.find("ATOM 1 C CA . ALA A 1");
    ASSERT_NE(atom_marker,std::string::npos);
    std::ofstream model_file(model_path);
    ASSERT_TRUE(model_file);
    model_file << template_text.substr(0,atom_marker);
    for(int index=0;index<25;++index)
    {
        const int serial=index+1;
        const int x=3+7*(index%5),y=3+7*(index/5);
        model_file << "ATOM " << serial << " C CA . ALA A " << serial << ' '
            << x << ' ' << y << " 3.0 1.0 0.0 1\n";
    }
    model_file << "ATOM 26 C CA . ALA B 1 3 6.5 3.0 1.0 0.0 1\n";
    model_file << "#\n";
    model_file.close();

    rgc::MapSimulationRequest simulation;
    simulation.model_file_path=model_path;
    simulation.output_dir=maps_dir;
    simulation.grid_spacing=.5;
    simulation.potential_model_choice=rgc::PotentialModel::SINGLE_GAUS;
    simulation.blurring_width_list={.55};
    simulation.verbosity=0;
    ASSERT_TRUE(rgc::RunCommand(simulation).succeeded);
    const auto map_path=FindGeneratedMap(maps_dir);
    ASSERT_FALSE(map_path.empty());

    rgc::PotentialAnalysisRequest request;
    request.model_file_path=model_path;
    request.map_file_path=map_path;
    request.selected_chain_ids={"A"};
    request.verbosity=0;

    struct FailureCase
    {
        const char * tag;
        rgc::detail::JointPostprocessingSnapshotHook hook;
        rgc::CommandOutcome outcome;
        int cli_exit;
    };
    for(const auto & failure:{
        FailureCase{"search-failure",&InjectJointSearchFailure,rgc::CommandOutcome::JointNotConverged,3},
        FailureCase{"no-donor",&InjectJointNoDonors,rgc::CommandOutcome::JointDownstreamUnavailable,4}})
    {
        SCOPED_TRACE(failure.tag);
        const bool qualified=failure.outcome==rgc::CommandOutcome::JointDownstreamUnavailable;
        const auto check_persistence=[&](const std::filesystem::path & database,const std::string & key) {
            SCOPED_TRACE(key);
            rg::DataRepository repository{database};
            const auto loaded=repository.LoadModel(key);
            ASSERT_NE(loaded,nullptr);
            const auto snapshot=loaded->GetAnalysisView().GetJointResult();
            ASSERT_TRUE(snapshot);
            ASSERT_TRUE(snapshot->initialization.valid);
            EXPECT_EQ(rg::IsQualifiedJointAnalysisResult(*snapshot),qualified);
            if(qualified)
            {
                EXPECT_TRUE(rg::HasUnavailableJointHaloMedianReconstruction(*snapshot));
                EXPECT_FALSE(snapshot->halo_median_reconstruction);
            }
            else
            {
                ASSERT_FALSE(snapshot->components.empty());
                EXPECT_EQ(snapshot->components.front().stop_reason,"block-recovery-exhausted");
            }
            for(const auto * atom:loaded->GetSelectedAtoms())
            {
                const auto view=rg::AtomLocalPotentialView::For(*atom);
                EXPECT_EQ(view.GetJointStageEstimate().point.has_value(),qualified);
                const auto & qscore=view.GetPostJointQScore();
                EXPECT_EQ(qscore.has_value(),qualified);
                if (qscore && qscore->missing_contributor_count)
                {
                    EXPECT_FALSE(qscore->scores);
                    EXPECT_EQ(qscore->status,rg::EvidenceStatus::Unavailable);
                    EXPECT_FALSE(qscore->reason.empty());
                }
                EXPECT_FALSE(view.GetPostFitPeeling());
                EXPECT_FALSE(view.GetPostJointLocalResult());
                EXPECT_FALSE(view.GetGroupMemberResult());
            }
        };

        request.saved_key_tag=failure.tag;
        request.database_path=directory.path()/(request.saved_key_tag+".sqlite");
        {
            ScopedJointPostprocessingSnapshotHook hook{failure.hook};
            const auto result=rgc::RunCommand(request);
            ASSERT_TRUE(result.succeeded);
            EXPECT_EQ(result.outcome,failure.outcome);
        }
        check_persistence(request.database_path,request.saved_key_tag);

        request.saved_key_tag+="-cli";
        request.database_path=directory.path()/(request.saved_key_tag+".sqlite");
        {
            ScopedJointPostprocessingSnapshotHook hook{failure.hook};
            EXPECT_EQ(RunCli({"RHBM-GEM","potential_analysis",
                "--model",model_path.string(),"--map",map_path.string(),
                "--database",request.database_path.string(),"--save-key",request.saved_key_tag,
                "--select-chain","A","--verbose","0"}),failure.cli_exit);
        }
        check_persistence(request.database_path,request.saved_key_tag);
    }
}

TEST(CommandApiPipelineTest, ProductionBoundedJointDefaultPersistsReloadsAndFeedsPostprocessing)
{
    command_test::ScopedTempDir directory{"bounded_joint_downstream"};
    const auto model_path=directory.path()/"bounded-input.cif";
    const auto template_path=command_test::TestDataPath("test_model.cif");
    std::ifstream template_file(template_path);
    ASSERT_TRUE(template_file);
    const std::string template_text{std::istreambuf_iterator<char>{template_file},
        std::istreambuf_iterator<char>{}};
    const auto atom_marker=template_text.find("ATOM 1 C CA . ALA A 1");
    ASSERT_NE(atom_marker,std::string::npos);
    std::ofstream model_file(model_path);
    ASSERT_TRUE(model_file);
    model_file << template_text.substr(0,atom_marker);
    for(int index=0;index<25;++index)
    {
        const int serial=index+1;
        const int x=index==0 ? 3 : index==1 ? 4 : 10+7*((index-2)%5);
        const int y=index<2 ? 3 : 3+7*((index-2)/5);
        model_file << "ATOM " << serial << " C CA . ALA A " << serial << ' '
            << x << ' ' << y << " 3.0 1.0 0.0 1\n";
    }
    model_file << "ATOM 26 C CA . ALA B 26 3.5 4.5 3.0 1.0 0.0 1\n";
    model_file << "#\n";
    model_file.close();

    const auto maps_dir=directory.path()/"maps";
    std::filesystem::create_directories(maps_dir);
    rgc::MapSimulationRequest simulation;
    simulation.model_file_path=model_path;
    simulation.output_dir=maps_dir;
    simulation.grid_spacing=.5;
    simulation.potential_model_choice=rgc::PotentialModel::SINGLE_GAUS;
    simulation.blurring_width_list={.55};
    simulation.verbosity=0;
    ASSERT_TRUE(rgc::RunCommand(simulation).succeeded);
    const auto map_path=FindGeneratedMap(maps_dir);
    ASSERT_FALSE(map_path.empty());

    const auto database_path=directory.path()/"bounded.sqlite";
    constexpr const char * key="bounded-joint";
    testing::internal::CaptureStdout();
    const auto cli_status=RunCli({"RHBM-GEM","potential_analysis",
        "--model",model_path.string(),"--map",map_path.string(),
        "--database",database_path.string(),"--save-key",key,
        "--select-chain","A","--verbose","4"});
    const auto joint_log=testing::internal::GetCapturedStdout();
    ASSERT_EQ(cli_status,0) << joint_log;
    EXPECT_NE(joint_log.find("Original-grid Peeling workspace="),std::string::npos);

    rg::DataRepository repository{database_path};
    const auto loaded=repository.LoadModel(key);
    ASSERT_NE(loaded,nullptr);
    const auto & stored=*loaded->GetAnalysisView().GetJointResult();
    auto reference_map=rg::ReadMap(map_path);
    auto reference_model=rg::ReadModel(model_path);
    ASSERT_NE(reference_map,nullptr);
    ASSERT_NE(reference_model,nullptr);
    ASSERT_TRUE(reference_map->MapValueArrayNormalization());
    reference_model->SelectAtoms([](const auto & atom) {return atom.GetChainID()=="A";});
    const auto reference=rgc::CaptureJointAnalysisResult(
        rgc::EstimateJointComponents(*reference_map,*reference_model));
    ASSERT_TRUE(rg::IsQualifiedJointAnalysisResult(reference));
    ExpectJointNumericalResultPreserved(reference,stored);
    ASSERT_TRUE(rg::IsQualifiedJointAnalysisResult(stored));
    ASSERT_TRUE(stored.selection_domain);
    EXPECT_EQ(stored.selection_domain->contributor_policy,"pure-halo-support-exclusion-v1");
    ASSERT_TRUE(stored.halo_median_reconstruction);
    const auto & halo_reconstruction=*stored.halo_median_reconstruction;
    EXPECT_EQ(halo_reconstruction.contract,rg::HaloMedianReconstructionContractVersion);
    EXPECT_EQ(halo_reconstruction.neighbor_count,rg::GeometryLocalTargetMedianNeighborCount);
    EXPECT_EQ(halo_reconstruction.donor_policy,rg::GeometryLocalTargetMedianDonorPolicy);
    EXPECT_EQ(halo_reconstruction.joint_run_id,
        rg::AtomLocalPotentialView::For(*loaded->FindAtomPtr(1))
            .GetJointStageEstimate().source.run_id);
    ASSERT_EQ(halo_reconstruction.halos.size(),1u);
    EXPECT_EQ(halo_reconstruction.halos.front().source,rg::HaloMedianSource::GeometryLocalTargetMedian);
    EXPECT_EQ(halo_reconstruction.halos.front().donor_atom_ids.size(),halo_reconstruction.halos.front().donor_count);
    EXPECT_LE(halo_reconstruction.halos.front().donor_count,rg::GeometryLocalTargetMedianNeighborCount);
    EXPECT_EQ(stored.selection_domain->target_indices.size(),25u);
    EXPECT_EQ(stored.atom_ids.size(),26u);
    const auto halo_id=std::find(stored.atom_ids.begin(),stored.atom_ids.end(),"26");
    ASSERT_NE(halo_id,stored.atom_ids.end());
    const auto halo_index=static_cast<std::size_t>(std::distance(stored.atom_ids.begin(),halo_id));
    EXPECT_EQ(std::find(stored.selection_domain->target_indices.begin(),
        stored.selection_domain->target_indices.end(),halo_index),
        stored.selection_domain->target_indices.end());
    ASSERT_NE(loaded->FindAtomPtr(1),nullptr);
    ASSERT_NE(loaded->FindAtomPtr(2),nullptr);
    ASSERT_NE(loaded->FindAtomPtr(26),nullptr);
    const auto target_one=loaded->FindAtomPtr(1)->GetPosition();
    const auto target_two=loaded->FindAtomPtr(2)->GetPosition();
    const auto halo_position=loaded->FindAtomPtr(26)->GetPosition();
    EXPECT_DOUBLE_EQ(std::hypot(target_one[0]-target_two[0],target_one[1]-target_two[1],
        target_one[2]-target_two[2]),1.0);
    EXPECT_LT(std::hypot(target_one[0]-halo_position[0],target_one[1]-halo_position[1],
        target_one[2]-halo_position[2]),2.5);
    ASSERT_TRUE(stored.metadata.solver);
    EXPECT_EQ(stored.metadata.solver->search_method,"gauss-seidel-bcd");
    EXPECT_EQ(stored.metadata.solver->bcd_local_search,"active-set-lm");
    const auto calibration=std::find_if(stored.evidence.begin(),stored.evidence.end(),[](const auto & item) {
        return item.name=="global-width-calibration";
    });
    ASSERT_NE(calibration,stored.evidence.end());
    EXPECT_EQ(calibration->status,rg::JointCheckStatus::Passed);
    EXPECT_NE(calibration->reason.find("global-initial-stage-p95-k4"),std::string::npos);
    ASSERT_TRUE(calibration->value && calibration->threshold);
    EXPECT_NEAR(*calibration->threshold,4.0*(*calibration->value),1e-12);
    const auto donor_coverage=std::find_if(stored.evidence.begin(),stored.evidence.end(),[](const auto & item) {
        return item.name=="qualified-initial-stage-donor-coverage";
    });
    ASSERT_NE(donor_coverage,stored.evidence.end());
    ASSERT_TRUE(donor_coverage->value && donor_coverage->threshold);
    EXPECT_GE(*donor_coverage->value,*donor_coverage->threshold);
    ASSERT_FALSE(stored.components.empty());
    EXPECT_TRUE(std::any_of(stored.components.begin(),stored.components.end(),[](const auto & component) {
        return std::any_of(component.evidence.begin(),component.evidence.end(),[](const auto & item) {
            return item.name=="active-set-lm-subproblems" && item.status==rg::JointCheckStatus::Passed;
        });
    }));

    std::size_t target_count{},peeling_count{},local_result_count{},local_available_count{},group_member_count{};
    for(const auto * atom:loaded->GetSelectedAtoms())
    {
        const auto local=rg::AtomLocalPotentialView::For(*atom);
        const auto & second=local.GetJointStageEstimate();
        if(second.source.method!=rg::EstimateMethod::JointComponents || second.source.role!=rg::FittingRole::Target)
            continue;
        ++target_count;
        EXPECT_TRUE(second.point);
        peeling_count+=local.GetPostFitPeeling().has_value();
        local_result_count+=local.GetPostJointLocalResult().has_value();
        if(local.GetPostJointLocalResult() && local.GetPostJointLocalResult()->status==rg::EvidenceStatus::Available)
            ++local_available_count;
        group_member_count+=local.GetGroupMemberResult().has_value();
    }
    EXPECT_EQ(target_count,25u);
    EXPECT_EQ(peeling_count,target_count);
    EXPECT_EQ(local_result_count,target_count);
    EXPECT_GT(local_available_count,0u);
    EXPECT_GT(group_member_count,0u);

    std::vector<std::pair<int, double>> valid_tail_ratios;
    std::size_t tail_raw_sample_count{};
    std::size_t tail_valid_sample_count{};
    std::size_t complete_tail_target_count{};
    for(const auto * atom:loaded->GetSelectedAtoms())
    {
        const auto view=rg::AtomLocalPotentialView::For(*atom);
        const auto assessment=view.AssessLocalFittingPeelingRatio(1.0,2.0);
        tail_raw_sample_count+=assessment.total_samples;
        tail_valid_sample_count+=assessment.valid_samples;
        if(assessment.total_samples && assessment.valid_samples==assessment.total_samples)
            ++complete_tail_target_count;
        if(assessment.ratio)
        {
            EXPECT_EQ(assessment.valid_samples,assessment.total_samples);
            valid_tail_ratios.emplace_back(atom->GetSerialID(),*assessment.ratio);
        }
    }
    ASSERT_GE(valid_tail_ratios.size(),3u);
    EXPECT_EQ(complete_tail_target_count,25u);
    EXPECT_EQ(tail_valid_sample_count,tail_raw_sample_count);
    EXPECT_EQ(valid_tail_ratios.size(),25u);
    RecordProperty("joint_target_count",static_cast<int>(target_count));
    RecordProperty("halo_contributor_count",static_cast<int>(stored.atom_ids.size()-target_count));
    RecordProperty("tail_complete_targets",static_cast<long long>(complete_tail_target_count));
    RecordProperty("tail_raw_samples",static_cast<long long>(tail_raw_sample_count));
    RecordProperty("tail_valid_samples",static_cast<long long>(tail_valid_sample_count));
    RecordProperty("tail_ratio_available_targets",static_cast<long long>(valid_tail_ratios.size()));

#ifdef RHBM_GEM_ENABLE_UMAP
    rgc::UmapEmbeddingRequest umap;
    umap.database_path=database_path;
    umap.model_key_tag=key;
    umap.output_dir=directory.path()/"umap";
    umap.num_neighbors=5;
    umap.num_epochs=5;
    umap.random_seed=42;
    umap.verbosity=0;
    ASSERT_TRUE(rgc::RunCommand(umap).succeeded);
    const auto csv_path=umap.output_dir/"umap_embedding_bounded-joint.csv";
    std::ifstream csv_file(csv_path);
    ASSERT_TRUE(csv_file);
    std::string csv_contents{std::istreambuf_iterator<char>{csv_file},std::istreambuf_iterator<char>{}};
    ASSERT_NE(csv_contents.find("tail peeling ratio"),std::string::npos);
    const auto first_newline=csv_contents.find('\n');
    ASSERT_NE(first_newline,std::string::npos);
    std::size_t embedded_row_count{};
    for(std::size_t offset=first_newline+1;offset<csv_contents.size();++offset)
        if(csv_contents[offset]=='\n') ++embedded_row_count;
    EXPECT_EQ(embedded_row_count,25u);
    RecordProperty("umap_included_rows",static_cast<long long>(embedded_row_count));
    const auto metadata_path=csv_path.string()+".metadata.json";
    std::ifstream metadata_file(metadata_path);
    ASSERT_TRUE(metadata_file);
    const std::string metadata_text{std::istreambuf_iterator<char>{metadata_file},std::istreambuf_iterator<char>{}};
    const auto metadata=boost::json::parse(metadata_text);
    EXPECT_EQ(metadata.at("estimator").as_string(),"joint-components");
    EXPECT_EQ(metadata.at("peeling_mode").as_string(),"original-grid-v1");
    EXPECT_EQ(metadata.at("included_rows").as_int64(),static_cast<std::int64_t>(embedded_row_count));

    ASSERT_TRUE(rgc::RunCommand(umap).succeeded);
    std::ifstream repeated_csv_file(csv_path);
    ASSERT_TRUE(repeated_csv_file);
    const std::string repeated_csv{std::istreambuf_iterator<char>{repeated_csv_file},std::istreambuf_iterator<char>{}};
    EXPECT_EQ(repeated_csv,csv_contents);
    RecordProperty("umap_csv_reproducible",repeated_csv==csv_contents ? "true" : "false");
#endif

    const auto reloaded_again=repository.LoadModel(key);
    ASSERT_NE(reloaded_again,nullptr);
    for(const auto & [serial_id,expected_ratio]:valid_tail_ratios)
    {
        const auto * atom=reloaded_again->FindAtomPtr(serial_id);
        ASSERT_NE(atom,nullptr);
        const auto ratio=rg::AtomLocalPotentialView::For(*atom).GetLocalFittingPeelingRatio(1.0,2.0);
        ASSERT_TRUE(ratio);
        EXPECT_DOUBLE_EQ(*ratio,expected_ratio);
    }

    const auto small_maps_dir=directory.path()/"small-maps";
    std::filesystem::create_directories(small_maps_dir);
    auto small_simulation=simulation;
    small_simulation.model_file_path=template_path;
    small_simulation.output_dir=small_maps_dir;
    ASSERT_TRUE(rgc::RunCommand(small_simulation).succeeded);
    const auto small_map_path=FindGeneratedMap(small_maps_dir);
    ASSERT_FALSE(small_map_path.empty());
    const auto small_database_path=directory.path()/"small-bounded.sqlite";
    constexpr const char * small_key="small-bounded-joint";
    int small_cli_status{};
    small_cli_status=RunCli({"RHBM-GEM","potential_analysis",
        "--model",template_path.string(),"--map",small_map_path.string(),
        "--database",small_database_path.string(),"--save-key",small_key,
        "--select-chain","A","--verbose","0"});
    EXPECT_NE(small_cli_status,0);
    rg::DataRepository small_repository{small_database_path};
    const auto small_loaded=small_repository.LoadModel(small_key);
    ASSERT_NE(small_loaded,nullptr);
    const auto & small_stored=*small_loaded->GetAnalysisView().GetJointResult();
    EXPECT_FALSE(rg::IsQualifiedJointAnalysisResult(small_stored));
    const auto unavailable_calibration=std::find_if(small_stored.evidence.begin(),small_stored.evidence.end(),
        [](const auto & item) {return item.name=="global-width-calibration";});
    ASSERT_NE(unavailable_calibration,small_stored.evidence.end());
    EXPECT_EQ(unavailable_calibration->status,rg::JointCheckStatus::Unavailable);
    for(const auto * atom:small_loaded->GetSelectedAtoms())
    {
        const auto local=rg::AtomLocalPotentialView::For(*atom);
        EXPECT_FALSE(local.GetJointStageEstimate().point);
        EXPECT_FALSE(local.GetPostFitPeeling());
        EXPECT_FALSE(local.GetPostJointLocalResult());
        EXPECT_FALSE(local.GetGroupMemberResult());
    }
}
