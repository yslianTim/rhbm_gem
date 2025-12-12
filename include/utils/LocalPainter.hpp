#pragma once

#include <string>
#include <vector>
#include <memory>
#include <Eigen/Dense>

class LocalPainter
{
public:
    LocalPainter(void);
    ~LocalPainter();

    static void PaintTemplate1(
        const Eigen::MatrixXd & data_matrix,
        const std::vector<double> & alpha_list,
        const std::string & x_axis_title,
        const std::string & file_name
    );

    static void PaintHistogram2D(
        const Eigen::MatrixXd & data_matrix,
        const std::string & x_axis_title,
        const std::string & y_axis_title,
        const std::string & file_name
    );

};
