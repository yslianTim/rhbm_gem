#include "FileFingerprint.hpp"
#include <boost/hash2/sha2.hpp>
#include <array>
#include <fstream>
#include <stdexcept>

namespace rhbm_gem {
std::string FileSha256(const std::filesystem::path & path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read file for SHA-256: " + path.string());
    boost::hash2::sha2_256 hash;
    std::array<char, 65536> buffer{};
    while (input.read(buffer.data(), static_cast<std::streamsize>(buffer.size())) || input.gcount() > 0)
        hash.update(buffer.data(), static_cast<size_t>(input.gcount()));
    if (!input.eof() || input.bad())
        throw std::runtime_error("Failed while hashing file: " + path.string());
    return boost::hash2::to_string(hash.result());
}

}
