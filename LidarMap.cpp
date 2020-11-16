#include "LidarMap.hpp"
#include <stdio.h>
#include <iostream>
#include <numeric>

const double LidarMap::MAP0[][2][2] = {
    {{   0.0 ,   0.0},{3000.0,    0.0}},
    {{   0.0,    0.0},{   0.0, 3000.0}},
    {{3000.0,    0.0},{3000.0, 3000.0}},
    {{   0.0, 3000.0},{3000.0, 3000.0}},
    {{1500.0, 1500.0},{3000.0, 1500.0}}
};

const double LidarMap::MAP2[][2][2] = {

};


LidarMap::LidarMap(unsigned int map_num){

    int segment_num = 0;
    switch (map_num){
    case 0:
        segment_num = std::extent<decltype(MAP0), 0>::value;
        setSegments(segment_num, MAP0);
        break;
    case 2:
        segment_num = std::extent<decltype(MAP2), 0>::value;
        setSegments(segment_num, MAP2);
        break;
    default:
        printf("[LidarMap] error map_num:%u\n", map_num);
        break;
    }
}

void LidarMap::setSegments(unsigned int segment_num, const double MAP[][2][2]){

    for (size_t i = 0; i < segment_num; i++){
        // mLineSegments.push_back(LineSegment(MAP[i][0][0], MAP[i][0][1], MAP[i][1][0], MAP[i][1][1]));
        mLineSegments.emplace_back(MAP[i][0][0], MAP[i][0][1], MAP[i][1][0], MAP[i][1][1]);
    }
}

std::vector<LineSegment> LidarMap::getLineSegments(){
    return mLineSegments;
}

std::vector<LineSegment> LidarMap::getLineSegments(const Eigen::Matrix2d& R, const Eigen::Vector2d& T){
    return calcLineSegments(mLineSegments, R, T);
}

std::vector<LineSegment> LidarMap::calcLineSegments(const std::vector<LineSegment>& lineSegments, const Eigen::Matrix2d& R, const Eigen::Vector2d& T){
    std::vector<LineSegment> movedLineSegments;
    movedLineSegments.reserve(lineSegments.size());

    for (const LineSegment &l : lineSegments){
        Eigen::Vector2d movedStart = R * l.start + T;
        Eigen::Vector2d movedEnd   = R * l.end   + T;
        // movedLineSegments.push_back(LineSegment(movedStart, movedEnd));
        movedLineSegments.emplace_back(movedStart, movedEnd);
    }
    
    return movedLineSegments;
}

std::pair<std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>, std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>> LidarMap::calcNearestPointsInMap2(const std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>& points, const std::vector<LineSegment>& lineSegments, double outlier_rate){
    
    size_t point_num = points.size();
    
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> nearest_points;
    nearest_points.reserve(point_num);

    std::vector<double> errors;
    errors.reserve(point_num);

    for (const Eigen::Vector2d &p : points){
        std::pair<Eigen::Vector2d, double> pair = calcNearestPointInMap2(p, lineSegments);
        nearest_points.push_back(pair.first);
        errors.push_back(pair.second);
    }
    
    double sum = std::accumulate(errors.begin(), errors.end(), 0.0);
    double mean = sum / errors.size();

    double sq_sum = std::inner_product(errors.begin(), errors.end(), errors.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / errors.size() - mean * mean);

    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> nearest_points_filtered;
    nearest_points_filtered.reserve(point_num);

    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> points_filtered;
    points_filtered.reserve(point_num);

    for (size_t i = 0; i < point_num; i++){
        if(errors[i] - mean < outlier_rate*stdev){
            nearest_points_filtered.push_back(nearest_points[i]);
            points_filtered.push_back(points[i]);
        }
    }

    return std::make_pair(nearest_points_filtered, points_filtered);
}

std::pair<Eigen::Vector2d, double> LidarMap::calcNearestPointInMap2(const Eigen::Vector2d& point, const std::vector<LineSegment>& lineSegments){
    std::vector<Eigen::Vector2d> near_points;
    near_points.reserve(lineSegments.size());

    std::vector<double> errors;
    errors.reserve(lineSegments.size());

    for (const LineSegment &l: lineSegments){
        Eigen::Vector2d near_point = calcNearestPoint(point, l);
        near_points.push_back(near_point);
        errors.push_back((point-near_point).norm());
    }

    size_t index = std::distance(errors.begin(), std::min_element(errors.begin(), errors.end()));
    return std::make_pair(near_points[index], errors[index]);
}

std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> LidarMap::calcNearestPointsInMap(const std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>& points, const std::vector<LineSegment>& lineSegments){
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> nearest_points;
    nearest_points.reserve(points.size());

    for (const Eigen::Vector2d &p : points){
        nearest_points.push_back(calcNearestPointInMap(p, lineSegments));
    }
    
    return nearest_points;
}

Eigen::Vector2d LidarMap::calcNearestPointInMap(const Eigen::Vector2d& point, const std::vector<LineSegment>& lineSegments){
    std::vector<Eigen::Vector2d> near_points;
    near_points.reserve(lineSegments.size());

    std::vector<double> errors;
    errors.reserve(lineSegments.size());

    for (const LineSegment &l: lineSegments){
        Eigen::Vector2d near_point = calcNearestPoint(point, l);
        near_points.push_back(near_point);
        errors.push_back((point-near_point).norm());
    }

    return near_points[std::distance(errors.begin(), std::min_element(errors.begin(), errors.end()))];
}

Eigen::Vector2d LidarMap::calcNearestPoint(const Eigen::Vector2d& point, const LineSegment& lineSegment){
    Eigen::Vector2d vecA = ls2Vec(lineSegment);
    Eigen::Vector2d vecB = pts2Vec(lineSegment.start, point);
    
    double dot_product = vecA.dot(vecB);

    if(dot_product > 0){
        double vecA_len = vecA.norm();
        double proj_len = dot_product / vecA_len;

        if(proj_len < vecA_len){
            Eigen::Vector2d proj_vec = vecA * (proj_len / vecA_len);
            return lineSegment.start + proj_vec;
        }else{
            return lineSegment.end;
        }
    }else{
        return lineSegment.start;
    }
}

Eigen::Vector2d LidarMap::ls2Vec(const LineSegment& lineSegment){
    return lineSegment.end - lineSegment.start; 
}

Eigen::Vector2d LidarMap::pts2Vec(const Eigen::Vector2d& start, const Eigen::Vector2d& end){
    return end - start;
}