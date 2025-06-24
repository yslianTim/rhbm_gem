#pragma once

#include <cstddef>
#include <tuple>
#include <vector>
#include <memory>
#include <unordered_map>
#include "GroupChargeEntry.hpp"

class AtomObject;
class ModelObject;
class AtomicChargeEntry;
enum class Element : uint16_t;
enum class Residue : uint16_t;
enum class Structure : uint8_t;

#ifdef HAVE_ROOT
class TH1D;
class TGraphErrors;
class TGraph2DErrors;
class TF1;
#endif

class ChargeEntryIterator
{
    AtomObject * m_atom_object;
    ModelObject * m_model_object;
    AtomicChargeEntry * m_atomic_entry;

public:
    ChargeEntryIterator(ModelObject * model_object);
    ChargeEntryIterator(AtomObject * atom_object);
    ~ChargeEntryIterator();
    double GetModelEstimateMinimum(int par_id, Element element) const;
    bool IsOutlierAtom(const std::string & class_key) const;
    bool IsAvailableGroupKey(uint64_t group_key, const std::string & class_key) const;
    size_t GetResidueCount(const std::string & class_key, Residue residue, Structure structure=static_cast<Structure>(0)) const;
    double GetModelEstimatePrior(uint64_t group_key, const std::string & class_key, int par_id) const;
    double GetModelVariancePrior(uint64_t group_key, const std::string & class_key, int par_id) const;
    const std::vector<AtomObject *> & GetAtomObjectList(uint64_t group_key, const std::string & class_key) const;
    std::unordered_map<int, AtomObject *> GetAtomObjectMap(uint64_t group_key, const std::string & class_key) const;

    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList(void) const;
    std::tuple<float, float> GetDistanceRange(double margin_rate=0.0) const;
    std::tuple<float, float> GetMapValueRange(double margin_rate=0.0) const;
    double GetModelEstimateMDPDE(int par_id) const;
    double GetInterceptEstimateMDPDE(void) const;
    double GetScaleEstimateMDPDE(void) const;
    double GetChargeEstimateMDPDE(void) const;

    #ifdef HAVE_ROOT
    std::unique_ptr<TH1D> CreateResidueCountHistogram(const std::string & class_key, Structure structure=static_cast<Structure>(0));
    std::unique_ptr<TGraphErrors> CreateModelEstimateToResidueGraph(std::vector<uint64_t> & group_key_list, const std::string & class_key, const int par_id=0);
    std::unique_ptr<TGraphErrors> CreateModelEstimateScatterGraph(std::vector<uint64_t> & group_key_list, const std::string & class_key, int par1_id=0, int par2_id=1);
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> CreateModelEstimateToResidueIDGraphMap(size_t main_chain_element_id, const int par_id=0, Residue residue=static_cast<Residue>(65535));
    std::unique_ptr<TGraphErrors> CreateModelBasisToMapValueGraph(int basis_id);
    std::unique_ptr<TGraphErrors> CreateRegressionDataGraph(int basis_id);
    #endif

private:
    bool IsAtomObjectAvailable(void) const;
    bool IsAtomicEntryAvailable(void) const;
    bool IsModelObjectAvailable(void) const;
    bool CheckGroupKey(uint64_t group_key, const std::string & class_key, bool verbose=true) const;
    Residue GetResidueFromGroupKey(uint64_t group_key, const std::string & class_key) const;

};
