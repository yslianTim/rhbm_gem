#include "core/command/detail/SimulationGeometry.hpp"
#include "core/estimator/joint_component/Problem.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <boost/json.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>

boost::json::array DiagnoseTruthReplay(const rhbm_gem::core::JointProblem &,
    const rhbm_gem::core::JointFitResult &,const rhbm_gem::core::joint_component::BcdExecutionPolicy &,bool);
void ExportBufferTruthGeometry(const rhbm_gem::MapObject &,const rhbm_gem::ModelObject &,
    const rhbm_gem::core::JointProblem &,const std::filesystem::path &);

namespace {
using namespace rhbm_gem;
namespace n=core::joint_component;
namespace j=boost::json;
namespace fs=std::filesystem;
using Point=std::array<double,3>;

double Number(const j::value & v)
{
    return v.to_number<double>();
}

Point ReadPoint(const j::value & value)
{
    const auto & a=value.as_array();
    return {Number(a.at(0)),Number(a.at(1)),Number(a.at(2))};
}

double Signal(double square,const Point & point)
{
    if (square>6.25) return 0.;
    const auto b=n::EvaluateKernel(square,point[1],2.5);
    return point[0]*b.gaussian+point[2]*b.charge;
}

struct Sample
{
    std::size_t row;
    std::vector<double> noise;
};

j::value SolveTruth(const j::object & input)
{
    const auto & name=input.at("name");
    std::vector<Point> positions,truth;
    for (const auto & p:input.at("positions").as_array()) positions.push_back(ReadPoint(p));
    for (const auto & p:input.at("truth").as_array()) truth.push_back(ReadPoint(p));
    if (positions.size()!=4 || truth.size()!=4) throw std::invalid_argument("Expected four fixture contributors.");
    const auto h=Number(input.at("spacing"));
    const Point spacing{h,h,h};
    const auto origin=ReadPoint(input.at("origin"));
    const auto & shape=input.at("dimensions").as_array();
    std::array<int,3> dimensions{};
    std::size_t count{1},map_bytes{};
    for (std::size_t i=0;i<dimensions.size();++i)
    {
        const auto dim=shape.at(i).as_int64();
        if (dim<=0 || dim>std::numeric_limits<int>::max() ||
            !n::CheckedMultiplyBytes(count,static_cast<std::size_t>(dim),count))
            throw std::invalid_argument("Invalid truth fixture dimensions.");
        dimensions[i]=static_cast<int>(dim);
    }
    if (!n::CheckedMultiplyBytes(count,2*sizeof(double),map_bytes) ||
        map_bytes>n::JointMemoryPlanningBudgetBytes)
        throw std::invalid_argument("Truth fixture exceeds the observation memory budget.");
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for (std::size_t i=0;i<positions.size();++i)
    {
        auto a=std::make_unique<AtomObject>();
        a->SetSerialID(static_cast<int>(i+1));
        a->SetElement(Element::CARBON);
        const auto & p=positions[i];
        a->SetPosition(p[0],p[1],p[2]);
        atoms.push_back(std::move(a));
    }
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const AtomObject & a){return a.GetSerialID()<=2;});
    std::vector<double> noiseless(count);
    for (int z=0;z<dimensions[2];++z) for (int y=0;y<dimensions[1];++y) for (int x=0;x<dimensions[0];++x)
    {
        const auto position=core::simulation::GridPosition({x,y,z},spacing,origin);
        const auto row=static_cast<std::size_t>(x)+static_cast<std::size_t>(dimensions[0])*
            (static_cast<std::size_t>(y)+static_cast<std::size_t>(dimensions[1])*static_cast<std::size_t>(z));
        for (std::size_t i=0;i<truth.size();++i)
            noiseless[row]+=Signal(core::simulation::SupportSquare(position,positions[i]),truth[i]);
    }
    std::vector<Sample> samples;
    for (const auto & v:input.at("noise").as_array())
    {
        const auto & s=v.as_array();
        const auto row=static_cast<std::size_t>(s.at(0).as_int64());
        if (row>=count) throw std::out_of_range("Noise row outside the fixture.");
        auto & entry=samples.emplace_back(Sample{row,{}});
        for (std::size_t i=1;i<s.size();++i) entry.noise.push_back(Number(s[i]));
    }
    const auto repeats=samples.at(0).noise.size();
    const auto * validation=input.if_contains("validation_noise");
    if (validation && validation->as_array().size()!=samples.size())
        throw std::invalid_argument("Independent reference noise must use the same voxel identities.");
    j::array results;
    std::vector<int> buffer{3};
    if (const auto * selection=input.if_contains("buffer"))
    {
        buffer.clear();
        for (const auto & id:selection->as_array()) buffer.push_back(static_cast<int>(id.as_int64()));
    }
    const bool width_diagnostics=input.if_contains("width_diagnostics");
    for (std::size_t rep=0;rep<repeats;++rep)
    {
        auto values=std::make_unique<double[]>(count);
        std::copy(noiseless.begin(),noiseless.end(),values.get());
        for (const auto & s:samples) values[s.row]+=s.noise.at(rep);
        const MapObject map(dimensions,spacing,origin,std::move(values));
        const auto baseline=core::BuildPureHaloSupportExclusionProblem(map,model);
        if (const auto * directory=input.if_contains("geometry_directory"))
        {
            ExportBufferTruthGeometry(map,model,baseline,fs::path(directory->as_string().c_str()));
            return j::object{{"scenario",name},{"geometry",*directory}};
        }
        const auto & common=baseline.FittingToOriginalRows();
        std::vector<std::size_t> recovered;
        for (std::string_view label:{"G0","single"})
        {
            if (label=="single" && buffer.empty())
            {
                auto identical=results.back().as_object();
                identical["strategy"]="single";
                identical["reused_identical_G0"]=true;
                results.push_back(std::move(identical));
                continue;
            }
            try
            {
                const auto start=std::chrono::steady_clock::now();
                const auto problem=label=="G0" ? baseline :
                    core::detail::BuildJointBufferProblem(map,model,buffer,
                        core::detail::BufferObservationDomain::Fixed).problem;
                const auto & fitting=problem.Input();
                if (label=="single")
                {
                    if (problem.OriginalInput().row_ids!=baseline.OriginalInput().row_ids ||
                        !std::includes(problem.FittingToOriginalRows().begin(),problem.FittingToOriginalRows().end(),
                            common.begin(),common.end()))
                        throw std::logic_error("Truth fixture violated fixed-domain row identity.");
                    std::set_difference(problem.FittingToOriginalRows().begin(),problem.FittingToOriginalRows().end(),
                        common.begin(),common.end(),std::back_inserter(recovered));
                }
                n::BcdExecutionPolicy policy;
                policy.width_bounds.lower=std::log(.25*h);
                policy.width_bounds.upper=std::log(Number(input.at("upper")));
                j::object recent_sweeps;
                n::JointProgressObserver observe;
                if (width_diagnostics)
                    observe=[&](const n::JointProgressEvent & event) {
                        if (!event.bcd) return;
                        const auto & b=*event.bcd;
                        auto & tail=recent_sweeps[event.component_id];
                        if (tail.is_null()) tail.emplace_array();
                        auto & a=tail.as_array();
                        if (a.size()==6) a.erase(a.begin());
                        a.push_back(j::object{{"sweep",b.sweep},{"objective",b.objective},
                            {"projected_width_kkt",b.projected_width_kkt},{"failed_blocks",b.failed_blocks},
                            {"recovery_attempts",b.recovery_attempts},{"recovery_successes",b.recovery_successes}});
                    };
                const auto fit=n::FitBcdComponents(problem,std::vector<double>(fitting.atom_ids.size(),.7),policy,observe);
                j::array checks;
                for (const auto & e:fit.evidence)
                    checks.push_back(j::object{{"name",e.name},{"status",std::string(joint_result_io::StatusText(e.status))},
                        {"reason",e.reason}});
                std::size_t sweeps{},evaluations{};
                j::array components;
                for (const auto & c:fit.components)
                {
                    sweeps+=c.search_sweeps;
                    evaluations+=static_cast<std::size_t>(c.profile_evaluations);
                    components.push_back(c.atoms.size());
                }
                const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                const bool passed=fit.RuntimeConvergence()==JointCheckStatus::Passed &&
                    fit.TargetRuntimeConvergence()==JointCheckStatus::Passed && fit.assembled_state.has_value();
                j::object row{{"replicate",rep},{"strategy",std::string(label)},{"qualified",passed},
                    {"checks",checks},{"evaluations",evaluations},{"sweeps",sweeps},{"component_sizes",components},
                    {"fitting_targets",fitting.atom_ids.size()},{"fitting_incidence",n::CountSupportEntries(fitting).value()},
                    {"construction_estimate",core::JointProblemAccess::Get(problem).construction_memory.estimated_peak_bytes},
                    {"wall_seconds",elapsed},{"recovered_rows",j::value_from(recovered)},
                    {"fitting_ids",j::value_from(fitting.atom_ids)}};
                if (width_diagnostics) row["recent_sweeps"]=std::move(recent_sweeps);
                if (width_diagnostics || (!passed && input.if_contains("diagnose")))
                {
                    try
                    {
                        row[width_diagnostics ? "width_diagnostics" : "uncertified_replay"]=
                            DiagnoseTruthReplay(problem,fit,policy,!width_diagnostics);
                    }
                    catch (const std::exception & error) {row["diagnostic_error"]=error.what();}
                }
                if (passed)
                {
                    const auto & state=*fit.assembled_state;
                    std::size_t zero_amplitude{},lower_face{},upper_face{};
                    std::vector<Point> coeff=truth;
                    j::array estimates;
                    for (std::size_t i=0;i<state.tau.size();++i)
                    {
                        const auto slot=static_cast<std::size_t>(std::stoi(fitting.atom_ids[i])-1);
                        coeff[slot]={state.ac[2*i],state.tau[i],state.ac[2*i+1]};
                        zero_amplitude+=state.ac[2*i]==0.0;
                        lower_face+=state.log_tau[i]-*policy.width_bounds.lower<=n::WidthBounds::RoundoffTolerance(*policy.width_bounds.lower);
                        upper_face+=*policy.width_bounds.upper-state.log_tau[i]<=n::WidthBounds::RoundoffTolerance(*policy.width_bounds.upper);
                        estimates.push_back(j::value_from(coeff[slot]));
                    }
                    row["bound_activity"]=j::array{zero_amplitude,lower_face,upper_face};
                    row["points"]=std::move(estimates);
                    const auto & original=problem.OriginalInput();
                    std::vector<double> prediction(original.row_ids.size());
                    // Simulated truth supplies the external Halo field only for the oracle prediction metric.
                    for (std::size_t i=0;i<original.atom_ids.size();++i)
                    {
                        const auto idx=static_cast<std::size_t>(std::stoi(original.atom_ids[i])-1);
                        for (const auto & sample:original.support[i])
                            prediction[sample.row]+=Signal(sample.squared_distance,coeff[idx]);
                    }
                    for (std::size_t i=0;i<prediction.size();++i)
                        prediction[i]-=noiseless[std::stoull(original.row_ids[i])];
                    auto score=[&](const auto & rows) -> j::value {
                        if (rows.empty()) return nullptr;
                        double sum{};
                        for (const auto i:rows) sum+=prediction[i]*prediction[i];
                        return std::sqrt(sum/static_cast<double>(rows.size()));
                    };
                    std::vector<std::size_t> all(original.row_ids.size()),core_rows;
                    std::iota(all.begin(),all.end(),0);
                    const auto first=std::find(original.atom_ids.begin(),original.atom_ids.end(),"1");
                    for (const auto & sample:original.support[static_cast<std::size_t>(first-original.atom_ids.begin())])
                        if (sample.squared_distance<=.25) core_rows.push_back(sample.row);
                    row["prediction_RMSE"]=j::object{{"original",score(all)},{"common",score(common)},
                        {"core",score(core_rows)}};
                    // Save per-row residuals only where either arm could report recovered RMSE.
                    j::array residuals;
                    for (const auto v:prediction) residuals.push_back(v);
                    row["residuals"]=std::move(residuals);
                    if (validation)
                    {
                        std::vector<double> separate_reference(count);
                        for (std::size_t i=0;i<samples.size();++i)
                        {
                            const auto & x=validation->as_array().at(i).as_array();
                            if (static_cast<std::size_t>(x.at(0).as_int64())!=samples[i].row)
                                throw std::invalid_argument("Independent reference noise row differs from training row.");
                            separate_reference[samples[i].row]=Number(x.at(rep+1));
                        }
                        j::array training_residual, validation_residual;
                        for (std::size_t i=0;i<prediction.size();++i)
                        {
                            const auto native_row=std::stoull(original.row_ids[i]);
                            training_residual.push_back(prediction[i]+noiseless[native_row]-map.GetMapValue(native_row));
                            validation_residual.push_back(prediction[i]-separate_reference[native_row]);
                        }
                        row["oracle_training_residuals"]=std::move(training_residual);
                        row["oracle_independent_residuals"]=std::move(validation_residual);
                        row["oracle_domain_rows"]=j::object{{"common",j::value_from(common)},
                            {"core",j::value_from(core_rows)}};
                    }
                }
                results.push_back(std::move(row));
            }
            catch (const std::exception & error)
            {
                results.push_back(j::object{{"replicate",rep},{"strategy",std::string(label)},
                    {"qualified",false},{"evaluations",0},{"error",error.what()}});
            }
        }
    }
    return j::object{{"scenario",name},{"truth",input.at("truth")},{"seed",input.at("seed")},
        {"spacing",h},{"upper",input.at("upper")},{"buffer",j::value_from(buffer)},{"repeats",repeats},{"fits",results}};
}
} // namespace

void RunBufferTruth(const std::filesystem::path & path, const std::filesystem::path & out)
{
    std::ifstream src(path);
    if (!src) throw std::runtime_error("Cannot read truth fixture.");
    const std::string payload{std::istreambuf_iterator<char>(src),{}};
    const auto input=boost::json::parse(payload);
    std::ofstream file(out);
    file.exceptions(std::ios::failbit | std::ios::badbit);
    for (const auto & scenario:input.as_array())
        file << boost::json::serialize(SolveTruth(scenario.as_object())) << '\n';
}
