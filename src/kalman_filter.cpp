#include "kalman_filter.h"
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::pow;
using std::sqrt;
using std::cout;
using std::endl;

#define PI 3.14159265

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {

  // measurement matrix - laser
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
             0, 1, 0, 0;

  // measurement matrix - jacobian
  H_jacobian_ = MatrixXd(3, 4);

  // state of object
  x_ = VectorXd(4);

  // state covariance matrix P
  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

  // state transition function
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0, // (0,2) to be replaced with delta_t at runtime
        0, 1, 0, 1, // (1,3) to be replaced with delta_t at runtime
        0, 0, 1, 0,
        0, 0, 0, 1;

}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - (H_laser_ * x_);
  MatrixXd S = (H_laser_ * P_ * H_laser_.transpose()) + R_;
  MatrixXd K = P_ * H_laser_.transpose() * S.inverse();
  
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_laser_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  cout << "Prediction px:" << px << ", py:" << py << ", vx:" << vx << ", vy:" << vy << endl;

  // check division by zero
  float c1 = px*px+py*py;
  if (fabs(c1) < 0.0001)
    return;
  // Check phi is not too close to pi or -pi
  if (fabs(PI - fabs(atan2(py, px))) < 0.005)
    return;

  VectorXd h_x(3);
  h_x << sqrt(px*px + py*py),
         atan2(py, px),
         (px*vx + py*vy)/sqrt(px*px + py*py);

  //cout << "atan2(py, px):" << atan2(py, px) << endl;

  VectorXd y = z - h_x;
  MatrixXd S = (H_jacobian_ * P_ * H_jacobian_.transpose()) + R_;
  MatrixXd K = P_ * H_jacobian_.transpose() * S.inverse();
  
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_jacobian_) * P_;

}
