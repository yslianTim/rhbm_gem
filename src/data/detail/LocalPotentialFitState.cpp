#include "data/detail/LocalPotentialFitState.hpp"

#include <stdexcept>
#include <utility>

namespace rhbm_gem {

LocalPotentialFitState::LocalPotentialFitState() = default;

LocalPotentialFitState::~LocalPotentialFitState() = default;

void LocalPotentialFitState::SetDataset(
    std::vector<Eigen::VectorXd> basis_and_response_entry_list)
{
    m_dataset = Dataset{ std::move(basis_and_response_entry_list) };
}

void LocalPotentialFitState::SetFitResult(FitResult value)
{
    m_fit_result = std::move(value);
}

const LocalPotentialFitState::Dataset & LocalPotentialFitState::GetDataset() const
{
    if (!m_dataset.has_value())
    {
        throw std::runtime_error("LocalPotentialFitState dataset is not available");
    }
    return *m_dataset;
}

const LocalPotentialFitState::FitResult & LocalPotentialFitState::GetFitResult() const
{
    if (!m_fit_result.has_value())
    {
        throw std::runtime_error("LocalPotentialFitState fit result is not available");
    }
    return *m_fit_result;
}

} // namespace rhbm_gem
