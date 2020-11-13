#include "ICP.hpp"

const double ICP::PI = 3.14159265358979323846;

ICP::ICP(double angle_range, double angle_offset):
    ANGLE_RANGE(angle_range),
    ANGLE_OFFSET(angle_offset)
{}

void ICP::estimate(const std::vector<double>& pts, LidarMap& lidarMap, const Eigen::Vector2d& initR, const Eigen::Vector2d& initT){
    std::vector<LineSegment> lineSegments = lidarMap.getLineSegments(initR, initT);
    // std::vector<LineSegment> lineSegments = lidarMap.getLineSegments();
}
 
std::vector<Eigen::Vector2d> ICP::dist_to_xy(const std::vector<double>& pts){
    int point_num = pts.size();

    std::vector<Eigen::Vector2d> pts_xy;
    pts_xy.reserve(point_num);

	for (int i = 0; i < point_num; i++) {
		double degree = (i / (double)(point_num - 1)) * ANGLE_RANGE + ANGLE_OFFSET;
		double theta = degToRad(degree);

        double x = pts[i] * std::cos(theta);
        double y = pts[i] * std::sin(theta);

        pts_xy.emplace_back(x, y);
	}

    return pts_xy;
}

std::vector<double> meter_to_mm(const std::vector<double>& pts){
    std::vector<double> new_pts(pts.size());
    std::copy(pts.begin(), pts.end(), new_pts.begin());

	std::for_each(new_pts.begin(), new_pts.end(), [](double& distance) { distance *= 1000.0; });    

    return new_pts;
}

double ICP::degToRad(double x) {
	return ((x)* PI / 180.0);
}
