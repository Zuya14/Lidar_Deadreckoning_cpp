#include <stdio.h>
#include "LidarMap.hpp"

#include <iostream>
#include <eigen3/Eigen/Eigen>

#include "ICP.hpp"

// #include <iostream>
// #include <stdio.h>
// #include <string>  
#include <fstream>  
// #include <sstream>   

#include <chrono>

void test() {
	std::string filename = "lidar_10_long_fixed.txt";
	// std::string filename = "lidar_25_long_fixed.txt";

	std::ifstream read_file(filename, std::ios::in);
	std::string read_line;

	if(!read_file){
		printf("can't find %s\n", filename.c_str());
		return;
	}	

	//printf("reading %s\n", filename.c_str());

	std::vector<double> points;
	points.reserve(1081);

	// LRF_Deadreckoning lrfDR;
    LidarMap lidarMap0 = LidarMap(0);
    ICP icp;

	double sum_time = 0.0;
	int count = 0;

    Eigen::Matrix2d R = Eigen::Matrix2d::Identity();
    Eigen::Vector2d T = Eigen::Vector2d::Zero();
    T << 1517, 752;

	while (!read_file.eof()) {
		points.clear();

		std::getline(read_file, read_line);
		std::string read_cell;
		std::istringstream iss(read_line);

		while (std::getline(iss, read_cell, '\t')) {
			points.push_back(std::stod(read_cell));
		}

		//printf("points.size:%u\n", points.size());

		if (points.size() == 1081) {
			std::chrono::system_clock::time_point start, end;
			start = std::chrono::system_clock::now();

			// bool success = lrfDR.update(points); 
			std::pair<Eigen::Matrix2d, Eigen::Vector2d> RT = icp.estimate(points, lidarMap0, R.transpose(), -T, 50, 1e-2);
            R = RT.first;
            T = RT.second;

            std::tuple<double, double, double> xyt = icp.calc_XY_Theta(R, T);
            double x =  std::get<0>(xyt);
            double y =  std::get<1>(xyt);
            double theta =  std::get<2>(xyt);
            
            end = std::chrono::system_clock::now();
			double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
            
            printf("(%d, %d, %1.2f), %f[ms]\n", (int)x, (int)y, theta, time);
            // printf("%f[ms]\n\n", time);
            sum_time += time;
            count++;
		}
	}

	printf("ave_time: %f[ms]\n\n", sum_time/count);
}

int main(){
    test();

    // LidarMap map0 = LidarMap(0);

    // Eigen::MatrixXd A(3, 2);
    // A << 1, 1,
    //     3, 3,
    //     5, 5;
    // Eigen::MatrixXd B(3, 2);
    // B << 10, 10,
    //     3, 3,
    //     1, 5;

    // std::cout << A << std::endl;
    // std::cout << B << std::endl;
    // std::cout << std::endl;
    // std::cout << A - B << std::endl;
    // std::cout << std::endl;

    // Eigen::VectorXd dist = (A - B).rowwise().norm();
    // double mean = dist.mean();
    // double std_dev = std::sqrt(dist.rowwise().squaredNorm().mean() - mean*mean);

    // std::cout << dist << std::endl;
    // std::cout << std::endl;
    // std::cout << mean << std::endl;
    // std::cout << std::endl;
    // std::cout << std_dev << std::endl;
    

    // double std_dev = std::sqrt((dist.).square().sum()/(dist.size()));

    return 0;
}