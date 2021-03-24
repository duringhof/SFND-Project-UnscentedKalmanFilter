#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

UKF::UKF() {

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd::Zero(n_x_);

  // create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);
  P_(3, 3) = std_laspx_ * std_laspy_;
  P_(4, 4) = std_laspx_ * std_laspy_;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  // Initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // Define spreading parameter
  lambda_ = 3 - n_aug_;

  // create vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  
  // set weights
  weights_(0) = lambda_ / (lambda_+n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; ++i) {

    weights_(i) = 1 / (2 * (lambda_ + n_aug_));
  }
}

UKF::~UKF() {}

// Generate augmented sigma points
void UKF::GenerateAugmentedSigmaPoints(MatrixXd &Xsig_aug) {
  
  // create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(7);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(7, 7);

  // create augmented mean state
  x_aug.head(n_x_) = x_;

  // create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
}

// Predict sigma points by using the process model
void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug, double delta_t) {

  // create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  // predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {

    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance() {

  // predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    x_diff(3) = Normalize(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

void UKF::PredictLidarMeasurement(VectorXd &z_pred, MatrixXd &S,
                                  MatrixXd &Zsig) {

  // transform sigma points into measurement space
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }

  // calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }

  // calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(2,2);
  R <<  std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
  S = S + R;
}

void UKF::PredictRadarMeasurement(VectorXd &z_pred, MatrixXd &S,
                                  MatrixXd &Zsig) {

  // transform sigma points into measurement space
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double rho = sqrt(p_x * p_x + p_y * p_y);
    double phi = atan2(p_y,p_x);
    double rhodot = (p_x * cos(yaw) * v + p_y * sin(yaw) * v) / rho;

    Zsig(0, i) = rho;
    Zsig(1, i) = phi;
    Zsig(2, i) = rhodot;
  }

  // calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }

  // calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { 
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    z_diff(1) = Normalize(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(3,3);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
  S = S + R;
}

void UKF::UpdateState(MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S,
                      VectorXd &z, MeasurementPackage::SensorType sensor_type) {

  // set measurement dimension
  int n_z;
  if (sensor_type == MeasurementPackage::RADAR) {
    // radar can measure r, phi, and r_dot
    n_z = 3;
  } else {
    // laser can measure px, py
    n_z = 2;
  }

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    z_diff(1) = Normalize(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    x_diff(3) = Normalize(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  z_diff(1) = Normalize(z_diff(1));

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) { // Initialize the state with the first measurement
    time_us_ = meas_package.timestamp_;

    // state vector x_: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and
    // rad
    double px, py, v = 0;
    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // set measurement dimension, radar can measure r, phi, and r_dot
      double r = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double r_dot = meas_package.raw_measurements_[2];
      px = r * cos(phi);
      py = r * sin(phi);
      if (abs(px) > 1e-6) {
        v = r_dot * sqrt(1 + pow((py / px), 2));
      }
      x_ << px, py, v, 0, 0;
      is_initialized_ = true;
    } else if (use_laser_ &&
               meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // set measurement dimension, laser can measure px, py
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];
      x_ << px, py, v, 0, 0;
      is_initialized_ = true;
    }
  } else { // Make prediction and update with new measurement
    double delta_t = (meas_package.timestamp_ - time_us_) / 1e6;

    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {

      double r = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double r_dot = meas_package.raw_measurements_[2];
      Prediction(delta_t);
      UpdateRadar(meas_package);
    } else if (use_laser_ &&
               meas_package.sensor_type_ == MeasurementPackage::LASER) {
      Prediction(delta_t);
      UpdateLidar(meas_package);
    }
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  GenerateAugmentedSigmaPoints(Xsig_aug);
  SigmaPointPrediction(Xsig_aug, delta_t);
  // Update state and covariance matrix by predicting it
  // using the process model
  PredictMeanAndCovariance();
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  time_us_ = meas_package.timestamp_;

  // Update state using lidar measurements
  // mean predicted measurement: lidar
  VectorXd z_pred_las = VectorXd::Zero(2);
  // measurement covariance matrix S
  MatrixXd S_las = MatrixXd::Zero(2, 2);
  // create matrix for sigma points in measurement space
  MatrixXd Zsig_las = MatrixXd(2, 2 * n_aug_ + 1);
  // predict lidar measurements
  PredictLidarMeasurement(z_pred_las, S_las, Zsig_las);

  UpdateState(Zsig_las, z_pred_las, S_las, meas_package.raw_measurements_,
              MeasurementPackage::LASER);
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  time_us_ = meas_package.timestamp_;
  // Update state using radar measurements
  // mean predicted measurement: radar
  VectorXd z_pred_rad = VectorXd::Zero(3);
  // measurement covariance matrix S
  MatrixXd S_rad = MatrixXd::Zero(3, 3);
  // create matrix for sigma points in measurement space
  MatrixXd Zsig_rad = MatrixXd(3, 2 * n_aug_ + 1);
  // predict radar measurements
  PredictRadarMeasurement(z_pred_rad, S_rad, Zsig_rad);

  UpdateState(Zsig_rad, z_pred_rad, S_rad, meas_package.raw_measurements_,
              MeasurementPackage::RADAR);
}

double UKF::Normalize(double angle) {

  while (angle > M_PI)
    angle -= 2. * M_PI;
  while (angle < -M_PI)
    angle += 2. * M_PI;
  return angle;
}