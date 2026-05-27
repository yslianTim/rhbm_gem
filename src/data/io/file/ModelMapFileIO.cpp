#include <rhbm_gem/data/io/ModelMapFileIO.hpp>

#include "CCP4Format.hpp"
#include "CifFormat.hpp"
#include "MrcFormat.hpp"
#include "PdbFormat.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

#include <array>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace rhbm_gem {

namespace {

using ModelReadFn = std::unique_ptr<ModelObject> (*)(std::istream &, const std::string &);
using ModelWriteFn = void (*)(std::ostream &, const ModelObject &, int);
using MapReadFn = std::unique_ptr<MapObject> (*)(std::istream &, const std::string &);
using MapWriteFn = void (*)(std::ostream &, const MapObject &);

struct ModelCodec
{
    std::string_view extension;
    ModelReadFn read;
    ModelWriteFn write;
};

struct MapCodec
{
    std::string_view extension;
    MapReadFn read;
    MapWriteFn write;
};

const auto & GetModelCodecs()
{
    static const std::array<ModelCodec, 4> codecs{
        ModelCodec{
            ".pdb",
            [](std::istream & file, const std::string & filename)
            {
                PdbFormat codec;
                return codec.ReadModel(file, filename);
            },
            [](std::ostream & outfile, const ModelObject & model_object, int model_parameter)
            {
                PdbFormat codec;
                codec.WriteModel(model_object, outfile, model_parameter);
            } },
        ModelCodec{
            ".cif",
            [](std::istream & file, const std::string & filename)
            {
                CifFormat codec;
                return codec.ReadModel(file, filename);
            },
            [](std::ostream & outfile, const ModelObject & model_object, int model_parameter)
            {
                CifFormat codec;
                codec.WriteModel(model_object, outfile, model_parameter);
            } },
        ModelCodec{
            ".mmcif",
            [](std::istream & file, const std::string & filename)
            {
                CifFormat codec;
                return codec.ReadModel(file, filename);
            },
            nullptr },
        ModelCodec{
            ".mcif",
            [](std::istream & file, const std::string & filename)
            {
                CifFormat codec;
                return codec.ReadModel(file, filename);
            },
            nullptr } };
    return codecs;
}

const auto & GetMapCodecs()
{
    static const std::array<MapCodec, 3> codecs{
        MapCodec{
            ".mrc",
            [](std::istream & file, const std::string & filename)
            {
                MrcFormat codec;
                return codec.ReadMap(file, filename);
            },
            [](std::ostream & outfile, const MapObject & map_object)
            {
                MrcFormat codec;
                codec.WriteMap(map_object, outfile);
            } },
        MapCodec{
            ".map",
            [](std::istream & file, const std::string & filename)
            {
                CCP4Format codec;
                return codec.ReadMap(file, filename);
            },
            [](std::ostream & outfile, const MapObject & map_object)
            {
                CCP4Format codec;
                codec.WriteMap(map_object, outfile);
            } },
        MapCodec{
            ".ccp4",
            [](std::istream & file, const std::string & filename)
            {
                CCP4Format codec;
                return codec.ReadMap(file, filename);
            },
            [](std::ostream & outfile, const MapObject & map_object)
            {
                CCP4Format codec;
                codec.WriteMap(map_object, outfile);
            } } };
    return codecs;
}

template <typename CodecList>
const typename CodecList::value_type & ResolveCodec(
    const CodecList & codecs,
    const std::filesystem::path & filename,
    bool for_write)
{
    const auto extension{ path_helper::GetExtension(filename) };
    for (const auto & codec : codecs)
    {
        if (codec.extension != extension)
        {
            continue;
        }
        if (for_write && codec.write == nullptr)
        {
            throw std::runtime_error("Unsupported file format for write: " + extension);
        }
        return codec;
    }
    throw std::runtime_error("Unsupported file format: " + extension);
}

template <typename StreamType>
StreamType OpenBinaryFile(const std::filesystem::path & filename, std::ios::openmode mode)
{
    StreamType stream{ filename, mode };
    if (!stream)
    {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    return stream;
}

} // namespace

std::unique_ptr<ModelObject> ReadModel(const std::filesystem::path & filename)
{
    try
    {
        const auto & codec{ ResolveCodec(GetModelCodecs(), filename, false) };
        auto file{ OpenBinaryFile<std::ifstream>(filename, std::ios::binary) };
        return codec.read(file, filename.string());
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to read model file '" + filename.string() + "': " + ex.what());
    }
}

void WriteModel(
    const std::filesystem::path & filename,
    const ModelObject & model_object,
    int model_parameter)
{
    try
    {
        const auto & codec{ ResolveCodec(GetModelCodecs(), filename, true) };
        auto outfile{ OpenBinaryFile<std::ofstream>(filename, std::ios::binary) };
        codec.write(outfile, model_object, model_parameter);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to write model file '" + filename.string() + "': " + ex.what());
    }
}

std::unique_ptr<MapObject> ReadMap(const std::filesystem::path & filename)
{
    try
    {
        const auto & codec{ ResolveCodec(GetMapCodecs(), filename, false) };
        auto file{ OpenBinaryFile<std::ifstream>(filename, std::ios::binary) };
        return codec.read(file, filename.string());
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to read map file '" + filename.string() + "': " + ex.what());
    }
}

void WriteMap(const std::filesystem::path & filename, const MapObject & map_object)
{
    try
    {
        const auto & codec{ ResolveCodec(GetMapCodecs(), filename, true) };
        auto outfile{ OpenBinaryFile<std::ofstream>(filename, std::ios::binary | std::ios::trunc) };
        codec.write(outfile, map_object);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to write map file '" + filename.string() + "': " + ex.what());
    }
}

} // namespace rhbm_gem
