#pragma once

#include <string>
#include <vector>
#include <memory>
#include <Eigen/Dense>

class LocalPainter
{
public:
    LocalPainter();
    ~LocalPainter();

    static void PaintTemplate1(
        const Eigen::MatrixXd & data_matrix,
        const std::vector<double> & alpha_list,
        const std::string & x_axis_title,
        const std::string & y_axis_title,
        const std::string & file_name
    );

    static void PaintHistogram2D(
        const std::vector<double> & data_array,
        int x_bin_size, double x_min, double x_max,
        int y_bin_size, double y_min, double y_max,
        const std::string & x_axis_title,
        const std::string & y_axis_title,
        const std::string & z_axis_title,
        const std::string & file_name
    );

};
