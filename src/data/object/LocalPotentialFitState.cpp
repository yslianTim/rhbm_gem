#include "data/object/LocalPotentialFitState.hpp"

#include <stdexcept>
#include <utility>

namespace rhbm_gem {

LocalPotentialFitState::LocalPotentialFitState() :
    m_sigma_square{ 0.0 }
{
}

LocalPotentialFitState::~LocalPotentialFitState() = default;

void LocalPotentialFitState::SetBasisAndResponseEntryList(
    std::vector<Eigen::VectorXd> value)
{
    m_basis_and_response_entry_list = std::move(value);
}

int LocalPotentialFitState::GetBasisAndResponseEntryListSize() const
{
    return static_cast<int>(m_basis_and_response_entry_list.size());
}

const std::vector<Eigen::VectorXd> & LocalPotentialFitState::GetBasisAndResponseEntryList() const
{
    if (m_basis_and_response_entry_list.empty())
    {
        throw std::runtime_error("LocalPotentialFitState basis/response list is empty");
    }
    return m_basis_and_response_entry_list;
}

} // namespace rhbm_gem
