#include "ICP.hpp"

const double ICP::PI = 3.14159265358979323846;

ICP::ICP(double angle_range, double angle_offset):
    ANGLE_RANGE(angle_range),
    ANGLE_OFFSET(angle_offset)
{}

std::pair<Eigen::Matrix2d, Eigen::Vector2d> ICP::estimate2(const std::vector<double>& pts, LidarMap& lidarMap, const Eigen::Matrix2d& initR, const Eigen::Vector2d& initT, size_t max_itr, double eps, double outlier_rate){
    std::vector<LineSegment> movedLS = lidarMap.getLineSegments(initR, initT);
    
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pts_xy = dist_to_xy(meter_to_mm(pts));

    double prev_mse = std::numeric_limits<double>::max();

    double max_outlier_rate = 3.5;
    if(outlier_rate > max_outlier_rate){
        outlier_rate = max_outlier_rate;
    }
	
    for (size_t iter = 0; iter < max_itr; iter++) {
        double rate = (outlier_rate - max_outlier_rate) * (iter/(double)max_itr) + max_outlier_rate;
        std::pair<std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>, std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>>> pts_pair = lidarMap.calcNearestPointsInMap2(pts_xy, movedLS, rate);
        
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> nearest_filtered = pts_pair.first;
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pts_xy_filtered = pts_pair.second;

        Eigen::MatrixXd mtx_pts_xy_filtered = convert_to_MatrixXd(pts_xy_filtered);
        Eigen::MatrixXd mtx_nearest_filtered = convert_to_MatrixXd(nearest_filtered);

        std::pair<Eigen::Matrix2d, Eigen::Vector2d> RT = fitTransform(mtx_nearest_filtered, mtx_pts_xy_filtered);
        Eigen::Matrix2d R = RT.first;
        Eigen::Vector2d T = RT.second;

        movedLS = lidarMap.calcLineSegments(movedLS, R, T);

        std::pair<double, double> error = calcError(mtx_nearest_filtered, mtx_pts_xy_filtered);
        double mse = error.first; 
        double std_dev = error.second;

        if(iter==0){
            // printf("T: %.1f, %.1f ", T[0], T[1]);
            // printf("%.1f, %.1f\n", mse, std_dev);
        }

        // printf("%lu\n", mtx_nearest_filtered.rows());

        if (std::abs(prev_mse - mse) < eps) {
            // printf("clear: %lu\n", iter);
            break;
        }

		prev_mse = mse;
    }

    std::vector<LineSegment> LS = lidarMap.getLineSegments();

    Eigen::MatrixXd mtx_LS = convert_to_MatrixXd(LS);
    Eigen::MatrixXd mtx_movedLS = convert_to_MatrixXd(movedLS);

    std::pair<Eigen::Matrix2d, Eigen::Vector2d> RT_final = fitTransform(mtx_movedLS, mtx_LS);

    return RT_final;
}

std::pair<Eigen::Matrix2d, Eigen::Vector2d> ICP::estimate(const std::vector<double>& pts, LidarMap& lidarMap, const Eigen::Matrix2d& initR, const Eigen::Vector2d& initT, size_t max_itr, double eps){
    std::vector<LineSegment> movedLS = lidarMap.getLineSegments(initR, initT);

    // std::vector<double> pts_mm = meter_to_mm(pts);
    // std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pts_xy = dist_to_xy(pts_mm);
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pts_xy = dist_to_xy(meter_to_mm(pts));
    Eigen::MatrixXd mtx_pts_xy = convert_to_MatrixXd(pts_xy);

    double prev_mse = std::numeric_limits<double>::max();

	for (size_t iter = 0; iter < max_itr; iter++) {
        std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> nearest = lidarMap.calcNearestPointsInMap(pts_xy, movedLS);

        Eigen::MatrixXd mtx_nearest = convert_to_MatrixXd(nearest);

        std::pair<Eigen::Matrix2d, Eigen::Vector2d> RT = fitTransform(mtx_nearest, mtx_pts_xy);
        Eigen::Matrix2d R = RT.first;
        Eigen::Vector2d T = RT.second;

        // if(iter==0){
        //     printf("T: %.1f, %.1f", T[0], T[1]);
        // }

        movedLS = lidarMap.calcLineSegments(movedLS, R, T);

        std::pair<double, double> error = calcError(mtx_nearest, mtx_pts_xy);
        double mse = error.first; 
        double std_dev = error.second;

        // printf("%f ", std::abs(prev_mse - mse));

        if (std::abs(prev_mse - mse) < eps) {
            // printf("clear: %lu\n", iter);
            break;
        }

		prev_mse = mse;
    }

    std::vector<LineSegment> LS = lidarMap.getLineSegments();

    Eigen::MatrixXd mtx_LS = convert_to_MatrixXd(LS);
    Eigen::MatrixXd mtx_movedLS = convert_to_MatrixXd(movedLS);

    std::pair<Eigen::Matrix2d, Eigen::Vector2d> RT_final = fitTransform(mtx_movedLS, mtx_LS);

    return RT_final;
}

std::tuple<double, double, double> ICP::calc_XY_Theta(const Eigen::Matrix2d& R, const Eigen::Vector2d& T){
    double x = T.x();
    double y = T.y();
    double thata = std::atan2(R(1, 0), R(0, 0));

    return std::make_tuple(x, y, thata);  
}

std::pair<double, double> ICP::calcError(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::VectorXd dist = (A - B).rowwise().norm();
    double mse = dist.mean();
    double std_dev = std::sqrt(dist.rowwise().squaredNorm().mean() - mse*mse);

    return std::make_pair(mse, std_dev);
}

std::pair<Eigen::Matrix2d, Eigen::Vector2d> ICP::fitTransform(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B){

	Eigen::Vector2d mean_A(0, 0), mean_B(0, 0);
	Eigen::MatrixXd AA = A, BB = B;

	mean_A = AA.colwise().mean();
	mean_B = BB.colwise().mean();

	AA.rowwise() -= mean_A.transpose();
	BB.rowwise() -= mean_B.transpose();

	Eigen::MatrixXd Sigma = AA.transpose() * BB;

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::MatrixXd U = svd.matrixU();
	Eigen::VectorXd S = svd.singularValues();
	Eigen::MatrixXd V = svd.matrixV();

	Eigen::Matrix2d R = V * U.transpose();

	Eigen::Vector2d T = mean_B - R * mean_A;    

    return std::make_pair(R, T);
}

Eigen::MatrixXd ICP::convert_to_MatrixXd(std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pts){

    size_t row = pts.size();
    Eigen::MatrixXd mtx = Eigen::MatrixXd(row, 2);

	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for (size_t i = 0; i < row; i++) {
		mtx(i, 0) = pts[i].x();
		mtx(i, 1) = pts[i].y();
	}

    return mtx;
}

Eigen::MatrixXd ICP::convert_to_MatrixXd(std::vector<LineSegment> lineSegments){

    size_t row = lineSegments.size();
    Eigen::MatrixXd mtx = Eigen::MatrixXd(row*2, 2);

    #ifdef _OPENMP
	#pragma omp parallel for
	#endif
	for (size_t i = 0; i < row; i++) {
        size_t idx1 = 2*i; 
        size_t idx2 = idx1 + 1;

		mtx(idx1, 0) = lineSegments[i].start.x();
		mtx(idx1, 1) = lineSegments[i].start.y();

        mtx(idx2, 0) = lineSegments[i].end.x();
		mtx(idx2, 1) = lineSegments[i].end.y();
	}

    return mtx;
}
 
std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> ICP::dist_to_xy(const std::vector<double>& pts){
    size_t point_num = pts.size();

    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> pts_xy;
    pts_xy.reserve(point_num);

	for (size_t i = 0; i < point_num; i++) {
		double degree = (i / (double)(point_num - 1)) * ANGLE_RANGE + ANGLE_OFFSET;
		double theta = degToRad(degree);

        double x = pts[i] * std::cos(theta);
        double y = pts[i] * std::sin(theta);

        pts_xy.emplace_back(x, y);
	}

    return pts_xy;
}

std::vector<double> ICP::meter_to_mm(const std::vector<double>& pts){
    std::vector<double> new_pts(pts.size());
    std::copy(pts.begin(), pts.end(), new_pts.begin());

	std::for_each(new_pts.begin(), new_pts.end(), [](double& distance) { distance *= 1000.0; });    

    return new_pts;
}

double ICP::degToRad(double x) {
	return ((x)* PI / 180.0);
}
