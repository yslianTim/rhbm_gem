#pragma once

#include <tuple>
#include <vector>
#include <memory>
#include <unordered_map>
#include "GroupPotentialEntry.hpp"

class AtomObject;
class ModelObject;
class AtomicPotentialEntry;
enum class Element : uint16_t;
enum class Residue : uint16_t;
enum class Structure : uint8_t;

#ifdef HAVE_ROOT
class TH1D;
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
    double GetGausEstimateMinimum(int par_id, Element element) const;
    bool IsOutlierAtom(const std::string & class_key) const;
    bool IsAvailableGroupKey(uint64_t group_key, const std::string & class_key) const;
    size_t GetResidueCount(const std::string & class_key, Residue residue, Structure structure=static_cast<Structure>(0)) const;
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
    std::unique_ptr<TH1D> CreateResidueCountHistogram(const std::string & class_key, Structure structure=static_cast<Structure>(0));
    std::vector<std::unique_ptr<TH1D>> CreateMainChainRankHistogram(int par_id, int & chain_size, Residue residue=static_cast<Residue>(65535), int extra_id=0, std::vector<Residue> veto_residues_list={});
    std::unique_ptr<TGraphErrors> CreateAmplitudeRatioToWidthScatterGraph(size_t target_id, size_t reference_id, Residue residue=static_cast<Residue>(65535));
    std::unique_ptr<TGraphErrors> CreateNormalizedGausEstimateScatterGraph(Element element, double reference_amplitude, bool reverse=false);
    std::unique_ptr<TGraphErrors> CreateBfactorToWidthScatterGraph(uint64_t group_key, const std::string & class_key);
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> CreateGausEstimateToResidueIDGraphMap(size_t main_chain_element_id, const int par_id=0, Residue residue=static_cast<Residue>(65535));
    std::unique_ptr<TGraphErrors> CreateGausEstimateToResidueGraph(std::vector<uint64_t> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateGausEstimateScatterGraph(std::vector<uint64_t> & group_key_list, const std::string & class_key, int par1_id=0, int par2_id=1);
    std::unique_ptr<TGraphErrors> CreateGausEstimateScatterGraph(Element element, bool reverse=false);
    std::unique_ptr<TGraphErrors> CreateGausEstimateScatterGraph(uint64_t group_key1, uint64_t group_key2, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateDistanceToMapValueGraph(void);
    std::unique_ptr<TGraphErrors> CreateBinnedDistanceToMapValueGraph(int bin_size=15, double x_min=0.0, double x_max=1.5);
    std::unique_ptr<TGraphErrors> CreateInRangeAtomsToGausEstimateGraph(uint64_t group_key, const std::string & class_key, double range=5.0, int par_id=0);
    std::unique_ptr<TGraphErrors> CreateCOMDistanceToGausEstimateGraph(uint64_t group_key, const std::string & class_key, int par_id=0);
    std::unique_ptr<TGraphErrors> CreateXYPositionTomographyGraph(double normalized_z_pos=0.5, double z_ratio_window=0.1, bool com_center=false);
    std::unique_ptr<TGraphErrors> CreateXYPositionTomographyGraph(std::vector<uint64_t> & group_key_list, const std::string & class_key, double normalized_z_pos=0.5, double z_ratio_window=0.1);
    std::unique_ptr<TGraph2DErrors> CreateXYPositionTomographyToGausEstimateGraph2D(std::vector<uint64_t> & group_key_list, const std::string & class_key, double normalized_z_pos=0.5, double z_ratio_window=0.1, int par_id=0);
    std::unique_ptr<TF1> CreateGroupGausFunctionPrior(uint64_t group_key, const std::string & class_key) const;
    #endif

private:
    bool IsAtomObjectAvailable(void) const;
    bool IsAtomicEntryAvailable(void) const;
    bool IsModelObjectAvailable(void) const;
    bool CheckGroupKey(uint64_t group_key, const std::string & class_key, bool verbose=true) const;
    Residue GetResidueFromGroupKey(uint64_t group_key, const std::string & class_key) const;

};