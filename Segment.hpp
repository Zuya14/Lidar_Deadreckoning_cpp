#ifndef SEGMENT_HPP_
#define SEGMENT_HPP_

#include "Point.hpp"

struct Segment{

    Segment(Point start, Point end):
        start(start), end(end)
    {}

    Segment(double x0, double y0, double x1, double y1){
        this(Point(x0, y0), Point(x1, y1));
    }

    Point start;
    Point end;
};

#endif /* SEGMENT_HPP_ */