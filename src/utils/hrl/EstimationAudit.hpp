#pragma once

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <chrono>
#include <cstddef>
#include <utility>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

namespace rhbm_gem::estimation_audit {

// Internal, observation-only context. Candidate workers inherit an owned scope.
struct Context
{
    std::size_t attempt{ 0 };
    std::string source;
    long atom{ -1 };
};
inline thread_local Context current;
inline bool Enabled() { return current.attempt >= 5 && current.attempt <= 8; }
class Scope
{
    Context previous{ current };
public:
    Scope(std::size_t attempt, std::string source, long atom = -1)
    {
#ifdef RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE
        current = { attempt, std::move(source), atom };
#else
        (void)attempt; (void)source; (void)atom;
#endif
    }
    ~Scope() { current = std::move(previous); }
};
inline std::string Number(double value)
{
    if (!std::isfinite(value)) return "null";
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}
inline void Emit(const std::string & kind, const std::string & fields) noexcept
{
    if (!Enabled()) return;
    try
    {
        std::ostringstream out;
        out << "Second-stage solver audit: schema=1, payload={\"attempt\":" << current.attempt
            << ",\"source\":" << std::quoted(current.source) << ",\"atom_index\":" << current.atom
            << ",\"kind\":" << std::quoted(kind) << ',' << fields << '}';
        Logger::Log(LogLevel::Debug, out.str());
    }
    catch (...) {} // Observation failures must not alter an estimator result.
}
} // namespace rhbm_gem::estimation_audit
