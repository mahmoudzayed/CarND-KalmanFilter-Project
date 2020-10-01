#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_;

  ///* Weights of sigma points
  double firstweight_;
  double remainingweights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

#ifdef DEBUG

  /// laser Chi-square 5% value for 2 dimension model
  double CHI_lidar;

  /// radar Chi-square 5% value for 2 dimension model
  double CHI_radar;
#endif

  ///The Normalized Innovation Squared for radar sensor
  std::vector<double> NIS_radar;

  ///The Normalized Innovation Squared for lidar sensor
  std::vector<double> NIS_lidar;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement 
   * Master module for controlling the follow of the filter.
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);  
  
  /**
   * initializeMeasurement
   * initialize the the filter states at the begining
   * @param meas_package The latest measurement data of either radar or laser
   */
  void initializeMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * augmanted matrix calculation
   * @param Xsig_aug augmanted matrix
   */
  void UKF::augmantaion(MatrixXd Xsig_aug);
  
  /**
   * predication model with a single
   * @param Xsig_aug sigma point
   */
  void Predictionmodel(VectorXd Xsig_aug, double delta_t);
		 
  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
#ifdef DEBUG
  /// Ch-square distribution 5% calculation va
  void ChiSquare_FIVEpec(MeasurementPackage meas_package, double CHI_value);
#endif
};

#endif /* UKF_H */
