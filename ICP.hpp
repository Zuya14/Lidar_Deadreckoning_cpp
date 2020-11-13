#ifndef ICP_HPP_
#define ICP_HPP_

#include <eigen3/Eigen/Eigen>
#include "LidarMap.hpp"

class ICP{
public:
    ICP(double angle_range=270.0, double angle_offset=-45.0);

    void estimate(const std::vector<double>& pts, LidarMap& lidarMap, const Eigen::Vector2d& initR, const Eigen::Vector2d& initT);

private:

    Eigen::Matrix2d fitTransform(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);

    std::vector<Eigen::Vector2d> dist_to_xy(const std::vector<double>& pts);
    std::vector<double> meter_to_mm(const std::vector<double>& pts);
    double degToRad(double x);

    const double ANGLE_RANGE;
    const double ANGLE_OFFSET;

    static const double PI;
};

#endif  /* ICP_HPP_ */