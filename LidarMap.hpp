#ifndef LIDAR_MAP_HPP_
#define LIDAR_MAP_HPP_

#include "LineSegment.hpp"
#include <vector>

class LidarMap{
public:
    LidarMap(unsigned int map_num);

    std::vector<LineSegment> getLineSegments();
    std::vector<LineSegment> getLineSegments(const Eigen::Matrix2d& R, const Eigen::Vector2d& T);

    std::vector<LineSegment> calcLineSegments(const std::vector<LineSegment>& lineSegments, const Eigen::Matrix2d& R, const Eigen::Vector2d& T);

    std::pair<std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector2d>> calcNearestPointsInMap2(const std::vector<Eigen::Vector2d>& points, const std::vector<LineSegment>& lineSegments);
    std::vector<Eigen::Vector2d> calcNearestPointsInMap(const std::vector<Eigen::Vector2d>& points, const std::vector<LineSegment>& lineSegments);

private:
    static const double MAP0[][2][2];
    static const double MAP2[][2][2];

    void setSegments(unsigned int segment_num, const double MAP[][2][2]);

    std::pair<Eigen::Vector2d, double> calcNearestPointInMap2(const Eigen::Vector2d& point, const std::vector<LineSegment>& lineSegments);
    Eigen::Vector2d calcNearestPointInMap(const Eigen::Vector2d& point, const std::vector<LineSegment>& lineSegments);
    Eigen::Vector2d calcNearestPoint(const Eigen::Vector2d& point, const LineSegment& lineSegment);

    Eigen::Vector2d ls2Vec(const LineSegment& lineSegment);
    Eigen::Vector2d pts2Vec(const Eigen::Vector2d& start, const Eigen::Vector2d& end);

    std::vector<LineSegment> mLineSegments;
};

#endif  /* LIDAR_MAP_HPP_ */