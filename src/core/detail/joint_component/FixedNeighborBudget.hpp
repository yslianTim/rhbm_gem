#pragma once
#include "ConvergenceCriteria.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace rhbm_gem::core::joint_component {
// Budget decisions never change the accepted state or qualify an endpoint.
class FixedNeighborBudget
{
    static constexpr std::size_t Window=20;
    std::array<double,Window> quality_{};
    std::size_t completed_{},limit_{},maximum_{};
    JointConvergenceCriteria criteria_;
    bool recovery_pending_{},extended_{};
public:
    FixedNeighborBudget(std::size_t nominal,std::size_t maximum,JointConvergenceCriteria criteria)
        :limit_(nominal ? std::min(nominal,maximum) : maximum),maximum_(maximum),criteria_(criteria) {}

    void Observe(double ac_kkt,double projected_width_kkt,bool recovery_pending)
    {
        double quality=std::numeric_limits<double>::infinity();
        if(IsValidJointConvergenceCriteria(criteria_) && std::isfinite(ac_kkt) && ac_kkt>=0 &&
            std::isfinite(projected_width_kkt) && projected_width_kkt>=0)
            quality=std::max(ac_kkt/criteria_.global_ac_kkt,
                projected_width_kkt/criteria_.global_width_gradient);
        quality_[completed_%Window]=quality;
        ++completed_;
        recovery_pending_=recovery_pending;
    }

    bool AllowsNextSweep()
    {
        if(completed_<limit_) return true;
        if(completed_>=maximum_ || completed_<Window || recovery_pending_) return false;
        double earlier=std::numeric_limits<double>::infinity(),recent=earlier;
        for(std::size_t i=0;i<Window;++i)
        {
            const double quality=quality_[(completed_-Window+i)%Window];
            if(!std::isfinite(quality)) return false;
            auto & best=i<Window/2 ? earlier : recent;
            best=std::min(best,quality);
        }
        // Best-of-half windows tolerate the measured rebounds and odd/even tail.
        // Only extend a nearby endpoint with at least 5% observed improvement.
        if(!(recent<=10.0 && recent<=.95*earlier && earlier>1.0)) return false;
        limit_=maximum_;
        extended_=true;
        return true;
    }

    bool HardLimitReached() const {return extended_ && completed_>=maximum_;}
};
}
