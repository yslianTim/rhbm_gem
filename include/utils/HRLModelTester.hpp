#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <Eigen/Dense>

class HRLModelTester
{
    const int m_basis_size;
    int m_member_size;
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> m_data_array;

public:
    HRLModelTester(void) = delete;
    explicit HRLModelTester(int basis_size, int member_size);
    ~HRLModelTester() = default;

    void BuildDataArray(void);

private:

};
