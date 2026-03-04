#include "FileFormatBackendFactory.hpp"

#include "CCP4Format.hpp"
#include "CifFormat.hpp"
#include "MapFileFormatBase.hpp"
#include "ModelFileFormatBase.hpp"
#include "MrcFormat.hpp"
#include "PdbFormat.hpp"

#include <stdexcept>

namespace rhbm_gem {

std::unique_ptr<ModelFileFormatBase> CreateModelFileFormatBackend(ModelFormatBackend backend)
{
    switch (backend)
    {
    case ModelFormatBackend::Pdb:
        return std::make_unique<PdbFormat>();
    case ModelFormatBackend::Cif:
        return std::make_unique<CifFormat>();
    }

    throw std::runtime_error("Unsupported model file backend.");
}

std::unique_ptr<MapFileFormatBase> CreateMapFileFormatBackend(MapFormatBackend backend)
{
    switch (backend)
    {
    case MapFormatBackend::Mrc:
        return std::make_unique<MrcFormat>();
    case MapFormatBackend::Ccp4:
        return std::make_unique<CCP4Format>();
    }

    throw std::runtime_error("Unsupported map file backend.");
}

} // namespace rhbm_gem
