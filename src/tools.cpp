#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::pow;
using std::sqrt;
using std::cout;
using std::endl;


Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << "Invalid input when calculating RMSE" << endl;
    return rmse;
  }

  for (int i=0; i < estimations.size(); ++i) {
    VectorXd diff = estimations[i] - ground_truth[i];
    VectorXd diff_2 = diff.array() * diff.array();
    rmse += diff_2;
  }

  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  // check division by zero
  float c1 = px*px+py*py;
  if (fabs(c1) < 0.0001) {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;
  }
  
  // compute the Jacobian matrix
  float px_2 = pow(px,2);
  float py_2 = pow(py,2);
  float vx_2 = pow(vx,2);
  float vy_2 = pow(vy,2);
  
  //Row 0
  Hj(0,0) = px/(sqrt(px_2 + py_2));
  Hj(0,1) = py/(sqrt(px_2 + py_2));
  Hj(0,2) = 0.0;
  Hj(0,3) = 0.0;

  //Row 1
  Hj(1,0) = -1.0 * (py/(px_2 + py_2));
  Hj(1,1) = px/(px_2 + py_2);
  Hj(1,2) = 0.0;
  Hj(1,3) = 0.0;

  //Row 2
  Hj(2,0) = py * ((vx*py - vy*px)/(sqrt(pow(px_2 + py_2,3))));
  Hj(2,1) = px * ((vy*px - vx*py)/(sqrt(pow(px_2 + py_2,3))));
  Hj(2,2) = px/(sqrt(px_2 + py_2));
  Hj(2,3) = py/(sqrt(px_2 + py_2));

  return Hj;
}
