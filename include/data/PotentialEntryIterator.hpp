#pragma once

#include <tuple>
#include <vector>
#include <unordered_map>
#include "GroupPotentialEntry.hpp"

class AtomObject;
class ModelObject;
class AtomicPotentialEntry;
enum class Element : uint16_t;

#ifdef HAVE_ROOT
class TGraphErrors;
class TGraph2DErrors;
class TF1;
#endif

class PotentialEntryIterator
{
    AtomObject * m_atom_object;
    ModelObject * m_model_object;
    AtomicPotentialEntry * m_atomic_entry;

public:
    PotentialEntryIterator(ModelObject * model_object);
    PotentialEntryIterator(AtomObject * atom_object);
    ~PotentialEntryIterator();
    bool IsAvailableGroupKey(uint64_t group_key, const std::string & class_key) const;
    std::tuple<double, double> GetGausEstimatePrior(uint64_t group_key, const std::string & class_key) const;
    std::tuple<double, double> GetGausVariancePrior(uint64_t group_key, const std::string & class_key) const;
    double GetGausEstimatePrior(uint64_t group_key, const std::string & class_key, int par_id) const;
    double GetGausVariancePrior(uint64_t group_key, const std::string & class_key, int par_id) const;
    const std::vector<AtomObject *> & GetAtomObjectList(uint64_t group_key, const std::string & class_key) const;
    std::unordered_map<int, AtomObject *> GetAtomObjectMap(uint64_t group_key, const std::string & class_key) const;

    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList(void) const;
    std::vector<std::tuple<float, float>> GetBinnedDistanceAndMapValueList(int bin_size=15, double x_min=0.0, double x_max=1.5) const;
    std::tuple<float, float> GetDistanceRange(double margin_rate=0.0) const;
    std::tuple<float, float> GetMapValueRange(double margin_rate=0.0) const;
    double GetAmplitudeEstimateMDPDE(void) const;
    double GetWidthEstimateMDPDE(void) const;

    #ifdef HAVE_ROOT
    std::unique_ptr<TGraphErrors> CreateGausEstimateToResidueGraph(std::vector<uint64_t> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateGausEstimateScatterGraph(std::vector<uint64_t> & group_key_list, const std::string & class_key, bool reverse=false);
    std::unique_ptr<TGraphErrors> CreateGausEstimateScatterGraph(Element element, bool reverse=false);
    std::unique_ptr<TGraphErrors> CreateGausEstimateScatterGraph(uint64_t group_key1, uint64_t group_key2, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateDistanceToMapValueGraph(void);
    std::unique_ptr<TGraphErrors> CreateBinnedDistanceToMapValueGraph(int bin_size=15, double x_min=0.0, double x_max=1.5);
    std::unique_ptr<TGraphErrors> CreateInRangeAtomsToGausEstimateGraph(uint64_t group_key, const std::string & class_key, double range=5.0, int par_id=0);
    std::unique_ptr<TGraphErrors> CreateXYPositionTomographyGraph(double normalized_z_pos=0.5, double z_ratio_window=0.1);
    std::unique_ptr<TGraphErrors> CreateXYPositionTomographyGraph(std::vector<uint64_t> & group_key_list, const std::string & class_key, double normalized_z_pos=0.5, double z_ratio_window=0.1);
    std::unique_ptr<TGraph2DErrors> CreateXYPositionTomographyToGausEstimateGraph2D(std::vector<uint64_t> & group_key_list, const std::string & class_key, double normalized_z_pos=0.5, double z_ratio_window=0.1, int par_id=0);
    std::unique_ptr<TF1> CreateGroupGausFunctionPrior(uint64_t group_key, const std::string & class_key) const;
    #endif

private:
    bool IsAtomObjectAvailable(void) const;
    bool IsAtomicEntryAvailable(void) const;
    bool IsModelObjectAvailable(void) const;
    bool CheckParameterIndex(int par_id) const;
    bool CheckGroupKey(uint64_t group_key, const std::string & class_key, bool verbose=true) const;

};