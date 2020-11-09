#ifndef ICP_HPP_
#define ICP_HPP_

#include <eigen3/Eigen/Eigen>

class ICP{
public:


private:

    Eigen::Matrix2d fitTransform(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);
};

#endif  /* ICP_HPP_ */