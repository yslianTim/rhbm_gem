#pragma once

#include <cstddef>
#include <string>
#include <vector>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>

#include <TAxis.h>

namespace rhbm_gem::painter_internal {

inline void RemodelAxisLabels(
    ::TAxis * axis,
    const std::vector<std::string> & label_list,
    double angle,
    int align)
{
    auto label_size{ static_cast<int>(label_list.size()) };
    root_helper::SetAxisTickAttribute(axis, 0.0, label_size + 1);
    axis->SetLimits(-0.75, static_cast<float>(label_size) - 0.5f);
    for (size_t i = 0; i < label_list.size(); i++)
    {
        auto label_index{ static_cast<int>(i) + 1 };
        axis->ChangeLabel(label_index, angle, -1, align, -1, -1, label_list.at(i).data());
    }
}

} // namespace rhbm_gem::painter_internal
#endif
