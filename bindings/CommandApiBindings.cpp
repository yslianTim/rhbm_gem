#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <rhbm_gem/core/command/CommandApi.hpp>

namespace py = pybind11;

namespace rhbm_gem::bindings {

void BindCommandApi(py::module_ & module)
{
    py::class_<CommonCommandRequest>(module, "CommonCommandRequest")
        .def(py::init<>())
        .def_readwrite("thread_size", &CommonCommandRequest::thread_size)
        .def_readwrite("verbose_level", &CommonCommandRequest::verbose_level)
        .def_readwrite("database_path", &CommonCommandRequest::database_path)
        .def_readwrite("folder_path", &CommonCommandRequest::folder_path);

    py::class_<PotentialAnalysisRequest>(module, "PotentialAnalysisRequest")
        .def(py::init<>())
        .def_readwrite("common", &PotentialAnalysisRequest::common)
        .def_readwrite("model_file_path", &PotentialAnalysisRequest::model_file_path)
        .def_readwrite("map_file_path", &PotentialAnalysisRequest::map_file_path)
        .def_readwrite("simulation_flag", &PotentialAnalysisRequest::simulation_flag)
        .def_readwrite("simulated_map_resolution", &PotentialAnalysisRequest::simulated_map_resolution)
        .def_readwrite("saved_key_tag", &PotentialAnalysisRequest::saved_key_tag)
        .def_readwrite("training_report_dir", &PotentialAnalysisRequest::training_report_dir)
        .def_readwrite("training_alpha_flag", &PotentialAnalysisRequest::training_alpha_flag)
        .def_readwrite("asymmetry_flag", &PotentialAnalysisRequest::asymmetry_flag)
        .def_readwrite("sampling_size", &PotentialAnalysisRequest::sampling_size)
        .def_readwrite("sampling_range_min", &PotentialAnalysisRequest::sampling_range_min)
        .def_readwrite("sampling_range_max", &PotentialAnalysisRequest::sampling_range_max)
        .def_readwrite("sampling_height", &PotentialAnalysisRequest::sampling_height)
        .def_readwrite("fit_range_min", &PotentialAnalysisRequest::fit_range_min)
        .def_readwrite("fit_range_max", &PotentialAnalysisRequest::fit_range_max)
        .def_readwrite("alpha_r", &PotentialAnalysisRequest::alpha_r)
        .def_readwrite("alpha_g", &PotentialAnalysisRequest::alpha_g);

    py::class_<PotentialDisplayRequest>(module, "PotentialDisplayRequest")
        .def(py::init<>())
        .def_readwrite("common", &PotentialDisplayRequest::common)
        .def_readwrite("painter_choice", &PotentialDisplayRequest::painter_choice)
        .def_readwrite("model_key_tag_list", &PotentialDisplayRequest::model_key_tag_list)
        .def_readwrite(
            "ref_model_key_tag_list",
            &PotentialDisplayRequest::ref_model_key_tag_list)
        .def_readwrite("pick_chain_id", &PotentialDisplayRequest::pick_chain_id)
        .def_readwrite("veto_chain_id", &PotentialDisplayRequest::veto_chain_id)
        .def_readwrite("pick_residue", &PotentialDisplayRequest::pick_residue)
        .def_readwrite("veto_residue", &PotentialDisplayRequest::veto_residue)
        .def_readwrite("pick_element", &PotentialDisplayRequest::pick_element)
        .def_readwrite("veto_element", &PotentialDisplayRequest::veto_element);

    py::class_<ResultDumpRequest>(module, "ResultDumpRequest")
        .def(py::init<>())
        .def_readwrite("common", &ResultDumpRequest::common)
        .def_readwrite("printer_choice", &ResultDumpRequest::printer_choice)
        .def_readwrite("model_key_tag_list", &ResultDumpRequest::model_key_tag_list)
        .def_readwrite("map_file_path", &ResultDumpRequest::map_file_path);

    py::class_<MapSimulationRequest>(module, "MapSimulationRequest")
        .def(py::init<>())
        .def_readwrite("common", &MapSimulationRequest::common)
        .def_readwrite("model_file_path", &MapSimulationRequest::model_file_path)
        .def_readwrite("map_file_name", &MapSimulationRequest::map_file_name)
        .def_readwrite(
            "potential_model_choice",
            &MapSimulationRequest::potential_model_choice)
        .def_readwrite(
            "partial_charge_choice",
            &MapSimulationRequest::partial_charge_choice)
        .def_readwrite("cutoff_distance", &MapSimulationRequest::cutoff_distance)
        .def_readwrite("grid_spacing", &MapSimulationRequest::grid_spacing)
        .def_readwrite("blurring_width_list", &MapSimulationRequest::blurring_width_list);

    py::class_<MapVisualizationRequest>(module, "MapVisualizationRequest")
        .def(py::init<>())
        .def_readwrite("common", &MapVisualizationRequest::common)
        .def_readwrite("model_file_path", &MapVisualizationRequest::model_file_path)
        .def_readwrite("map_file_path", &MapVisualizationRequest::map_file_path)
        .def_readwrite("atom_serial_id", &MapVisualizationRequest::atom_serial_id)
        .def_readwrite("sampling_size", &MapVisualizationRequest::sampling_size)
        .def_readwrite("window_size", &MapVisualizationRequest::window_size);

    py::class_<PositionEstimationRequest>(module, "PositionEstimationRequest")
        .def(py::init<>())
        .def_readwrite("common", &PositionEstimationRequest::common)
        .def_readwrite("map_file_path", &PositionEstimationRequest::map_file_path)
        .def_readwrite("iteration_count", &PositionEstimationRequest::iteration_count)
        .def_readwrite("knn_size", &PositionEstimationRequest::knn_size)
        .def_readwrite("alpha", &PositionEstimationRequest::alpha)
        .def_readwrite("threshold_ratio", &PositionEstimationRequest::threshold_ratio)
        .def_readwrite("dedup_tolerance", &PositionEstimationRequest::dedup_tolerance);

    py::class_<HRLModelTestRequest>(module, "HRLModelTestRequest")
        .def(py::init<>())
        .def_readwrite("common", &HRLModelTestRequest::common)
        .def_readwrite("tester_choice", &HRLModelTestRequest::tester_choice)
        .def_readwrite("fit_range_min", &HRLModelTestRequest::fit_range_min)
        .def_readwrite("fit_range_max", &HRLModelTestRequest::fit_range_max)
        .def_readwrite("alpha_r", &HRLModelTestRequest::alpha_r)
        .def_readwrite("alpha_g", &HRLModelTestRequest::alpha_g);

    py::class_<ExecutionReport>(module, "ExecutionReport")
        .def(py::init<>())
        .def_readonly("prepared", &ExecutionReport::prepared)
        .def_readonly("executed", &ExecutionReport::executed)
        .def_readonly("validation_issues", &ExecutionReport::validation_issues);

// BEGIN GENERATED: command-pybind-run-bindings
    module.def("RunPotentialAnalysis", &RunPotentialAnalysis);
    module.def("RunPotentialDisplay", &RunPotentialDisplay);
    module.def("RunResultDump", &RunResultDump);
    module.def("RunMapSimulation", &RunMapSimulation);
    module.def("RunMapVisualization", &RunMapVisualization);
    module.def("RunPositionEstimation", &RunPositionEstimation);
    module.def("RunHRLModelTest", &RunHRLModelTest);
// END GENERATED: command-pybind-run-bindings
}

} // namespace rhbm_gem::bindings
