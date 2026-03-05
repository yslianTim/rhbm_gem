#pragma once

#include "DataObjectVisitor.hpp"

namespace rhbm_gem {

struct PotentialAnalysisCommandOptions;
class AtomObject;
class BondObject;
class MapObject;
class ModelObject;

class AtomSamplingVisitor : public DataObjectVisitor
{
    const MapObject * m_map_object;
    const PotentialAnalysisCommandOptions * m_options;

public:
    AtomSamplingVisitor(
        const MapObject & map_object,
        const PotentialAnalysisCommandOptions & options);

    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;
};

class AtomGroupingVisitor : public DataObjectVisitor
{
public:
    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;
};

class LocalFittingVisitor : public DataObjectVisitor
{
    const PotentialAnalysisCommandOptions * m_options;
    double m_universal_alpha_r;

public:
    LocalFittingVisitor(
        const PotentialAnalysisCommandOptions & options,
        double universal_alpha_r);

    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;
};

class BondSamplingVisitor : public DataObjectVisitor
{
    const MapObject * m_map_object;
    const PotentialAnalysisCommandOptions * m_options;

public:
    BondSamplingVisitor(
        const MapObject & map_object,
        const PotentialAnalysisCommandOptions & options);

    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;
};

class BondGroupingVisitor : public DataObjectVisitor
{
public:
    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;
};

} // namespace rhbm_gem
