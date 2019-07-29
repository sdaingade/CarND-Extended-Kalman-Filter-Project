#ifndef FusionEKF_H_
#define FusionEKF_H_

#include <fstream>
#include <string>
#include <vector>
#include "Eigen/Dense"
#include "kalman_filter.h"
#include "measurement_package.h"
#include "tools.h"

using Eigen::MatrixXd;

class FusionEKF {
 public:
  /**
   * Constructor.
   */
  FusionEKF();

  /**
   * Destructor.
   */
  virtual ~FusionEKF();

  /**
   * Run the whole flow of the Kalman Filter from here.
   */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  MatrixXd SetupQ(float dt, float noise_ax, float noise_ay);

  /**
   * Kalman Filter update and prediction math lives in here.
   */
  KalmanFilter ekf_;

 private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;

  // measurement covariance matrix
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;

  // Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  float noise_ax = 9;
  float noise_ay = 9;
  
  int iter_ = 1;
};

#endif // FusionEKF_H_
