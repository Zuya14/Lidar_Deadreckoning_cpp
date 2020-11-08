#ifndef LINESEGMENT_HPP_
#define LINESEGMENT_HPP_

#include <eigen3/Eigen/Eigen>

struct LineSegment{

    LineSegment(Eigen::Vector2d start, Eigen::Vector2d end):
        start(start), end(end)
    {}

    LineSegment(double x0, double y0, double x1, double y1):
        LineSegment(Eigen::Vector2d(x0, y0), Eigen::Vector2d(x1, y1))
    {}

    Eigen::Vector2d start;
    Eigen::Vector2d end;
};

#endif /* LINESEGMENT_HPP_ */