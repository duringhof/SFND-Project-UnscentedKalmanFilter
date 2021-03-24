#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
 
  UKF();

  virtual ~UKF();

  void ProcessMeasurement(MeasurementPackage meas_package);
  void Prediction(double delta_t);
  void UpdateLidar(MeasurementPackage meas_package);
  void UpdateRadar(MeasurementPackage meas_package);

  void GenerateAugmentedSigmaPoints(MatrixXd &Xsig_aug);
  void SigmaPointPrediction(MatrixXd &Xsig_aug, double delta_t);
  void PredictMeanAndCovariance();
  void PredictRadarMeasurement(VectorXd &z_pred, MatrixXd &S, MatrixXd &Zsig);
  void PredictLidarMeasurement(VectorXd &z_pred, MatrixXd &S, MatrixXd &Zsig);
  void UpdateState(MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S, VectorXd &z,
                   MeasurementPackage::SensorType sensor_type);

  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  // state covariance matrix
  MatrixXd P_;

  // predicted sigma points matrix
  MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_;

  // Weights of sigma points
  VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Sigma point spreading parameter
  double lambda_;

  // radar measurement dimension
  int nrad_z_;

  // laser measurement dimension
  int nlas_z_;
};

#endif // UKF_H