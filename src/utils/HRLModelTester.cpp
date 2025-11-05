#include "HRLModelTester.hpp"
#include "HRLModelHelper.hpp"
#include "Logger.hpp"

#include <cmath>
#include <random>
#include <limits>
#include <utility>
#include <stdexcept>
#include <sstream>

namespace
{
int ValidatePositive(int value, const char * name)
{
    if (value <= 0)
    {
        throw std::invalid_argument(std::string(name) + " must be positive value");
    }
    return value;
}
} // namespace

HRLModelTester::HRLModelTester(int basis_size, int member_size) :
    m_basis_size{ ValidatePositive(basis_size, "basis_size") },
    m_member_size{ ValidatePositive(member_size, "member_size") },
    m_data_array{}
{
    m_data_array.reserve(static_cast<size_t>(m_member_size));
}

void HRLModelTester::BuildDataArray(void)
{
    m_data_array.clear();
}
