#ifndef ICP_HPP_
#define ICP_HPP_

#include <eigen3/Eigen/Eigen>
#include "LidarMap.hpp"

class ICP{
public:
    ICP(double angle_range=270.0, double angle_offset=-45.0);

    std::pair<Eigen::Matrix2d, Eigen::Vector2d> estimate(const std::vector<double>& pts, LidarMap& lidarMap, const Eigen::Matrix2d& initR, const Eigen::Vector2d& initT, size_t max_itr=10, double eps=1e-6);

    std::tuple<double, double, double> calc_XY_Theta(const Eigen::Matrix2d& R, const Eigen::Vector2d& T);

private:

    std::pair<double, double> calcError(Eigen::MatrixXd A, Eigen::MatrixXd B);

    std::pair<Eigen::Matrix2d, Eigen::Vector2d> fitTransform(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);

    Eigen::MatrixXd convert_to_MatrixXd(std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pts);
    Eigen::MatrixXd convert_to_MatrixXd(std::vector<LineSegment> lineSegments);

    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> dist_to_xy(const std::vector<double>& pts);
    std::vector<double> meter_to_mm(const std::vector<double>& pts);
    double degToRad(double x);

    const double ANGLE_RANGE;
    const double ANGLE_OFFSET;

    static const double PI;
};

#endif  /* ICP_HPP_ */