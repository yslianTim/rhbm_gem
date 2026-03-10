#pragma once

// Centralized seam for data-layer internal test-only dependencies.
#include "internal/io/file/FileFormatRegistry.hpp"
#include "internal/io/file/FileProcessFactoryBase.hpp"
#include "internal/io/file/FileProcessFactoryResolver.hpp"
#include "internal/io/sqlite/DataObjectDAOBase.hpp"
#include "internal/io/sqlite/DataObjectDAOFactoryRegistry.hpp"
#include "internal/io/sqlite/DatabaseManager.hpp"
#include "internal/io/sqlite/MapObjectDAO.hpp"
#include "internal/io/sqlite/ModelObjectDAO.hpp"
#include "internal/io/sqlite/SQLiteWrapper.hpp"
