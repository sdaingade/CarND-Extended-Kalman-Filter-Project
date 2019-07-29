#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <cmath>

#define PI 3.14159265

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  cout << "Iteration: " << iter_ << endl;
  iter_++;

  /**
   * Initialization
   */
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    // State covariance matrix is initialized in KalmanFilter::KalmanFilter()
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "Initialization using radar measurement" << endl;
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      
      // Calculate angle with y axis and resolve rho into x and y components
      float angle_with_y = (PI/2.0) - phi;
      ekf_.x_ << rho * sin(angle_with_y),
                 -1 * rho * cos(angle_with_y),
                 0,
                 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "Initialization using laser measurment" << endl;
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                 measurement_pack.raw_measurements_[1],
                 0,
                 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  // Update the state transition matrix F according to the new elapsed time.
  // Time is measured in seconds.
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Update the state transition matrix F according to the new elapsed time.
  // Time is measured in seconds.
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  // Update the process noise covariance matrix.
  // Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  ekf_.Q_ = SetupQ(dt, noise_ax, noise_ay);

  //cout << "Before calling predict" << endl;
  ekf_.Predict();
  //cout << "After calling predict" << endl;

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.R_ = R_radar_;
    //cout << "Before calling CalculateJacobian" << measurement_pack.raw_measurements_ << endl;
    ekf_.H_jacobian_ = tools.CalculateJacobian(ekf_.x_);
    //cout << "After calling CalculateJacobian" << endl;
    //cout << "Before calling UpdateEKF - RADAR" << endl;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    //cout << "After calling UpdateEKF - RADAR" << endl;
  } else {
    //cout << "Before calling Update - Lidar" << endl;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
    //cout << "After calling Update - Lidar" << endl;
  }

  // print the output
  cout << "After update x_ = " << ekf_.x_ << endl;
  cout << "After update P_ = " << ekf_.P_ << endl;
}

MatrixXd FusionEKF::SetupQ(float dt, float noise_ax, float noise_ay) {
    MatrixXd Q(4,4);

    //Row 0
    Q(0,0) = (pow(dt, 4)/4.0) * noise_ax;
    Q(0,1) = 0.0;
    Q(0,2) = (pow(dt, 3)/2.0) * noise_ax;
    Q(0,3) = 0.0;

    //Row 1
    Q(1,0) = 0.0;
    Q(1,1) = (pow(dt, 4)/4.0) * noise_ay;
    Q(1,2) = 0.0;
    Q(1,3) = (pow(dt, 3)/2.0) * noise_ay;

    //Row 2
    Q(2,0) = (pow(dt, 3)/2.0) * noise_ax;
    Q(2,1) = 0.0;
    Q(2,2) = pow(dt, 2) * noise_ax;
    Q(2,3) = 0.0;

    //Row 3
    Q(3,0) = 0.0;
    Q(3,1) = (pow(dt, 3)/2.0) * noise_ay;
    Q(3,2) = 0.0;
    Q(3,3) = pow(dt, 2) * noise_ay;

    return Q;    
}