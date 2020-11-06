#ifndef LIDAR_MAP_HPP_
#define LIDAR_MAP_HPP_

#include "Segment.hpp"

class LidarMap{
public:
    LidarMap(int map_num);

private:
    Segment mSegments;
};

#endif  /* LIDAR_MAP_HPP_ */