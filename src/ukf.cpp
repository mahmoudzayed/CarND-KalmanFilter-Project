#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 3;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.8;

	//DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;
	//DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

	// initialize the state vector  with the first measurement reading.
	is_initialized_ = false;

	///* State dimension
	n_x_ = 5;

	///* Augmented state dimension
	n_aug_ = 7;

	///* Sigma point spreading parameter
	lambda_ = 3 - n_aug_;

	///create matrix with predicted sigma points as columns
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	/// set weights
	firstweight_ = lambda_ / (lambda_ + n_aug_);
	remainingweights_ = 0.5 / (n_aug_ + lambda_);

	P_ = MatrixXd(n_x_, n_x_);

	x_ = VectorXd(n_x_);

#ifdef DEBUG
	CHI_radar = 7.815;

	CHI_lidar = 5.991;
#endif
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::initializeMeasurement(MeasurementPackage meas_package) {
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		/**
		Convert radar from polar to cartesian coordinates and initialize state.
		*/
		double range = meas_package.raw_measurements_[0];
		double bearing = meas_package.raw_measurements_[1];

		x_ << range * cos(bearing), range*sin(bearing), 0, 0, 0;

	}
	else
	{
		if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			//set the state with the initial location and zero velocity
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

		}
	}
	P_ = 0.1 * MatrixXd::Identity(n_x_, n_x_);
	time_us_ = meas_package.timestamp_;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
	Make sure you switch between lidar and radar measurements.
	*/
	if (!is_initialized_)
	{
		initializeMeasurement(meas_package);
		is_initialized_ = true;
		return;
	}

	/// calculate the time step for this measurment
	double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;

	/// predict the State vector using UKF 
	Prediction(delta_t);

	/// update the measurment and state vector
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		UpdateRadar(meas_package);
	}
	else
	{
		if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			UpdateLidar(meas_package);
		}
	}

}

/**
 * create augmantation matrix.
 * @param {MatrixXd} Xsig_aug the newly created augmantation matrix
 */
void UKF::augmantaion(MatrixXd Xsig_aug) {
	
	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);
	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	// calculate once to reduce calculation cost.
	double sqrtlambda = sqrt(lambda_ + n_aug_);
	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.fill(0.0);
	Xsig_aug.col(0) = x_aug;

	for (int i = 0; i < n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrtlambda * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrtlambda * L.col(i);
	}
}

/**
 * Predicts model for sigma points
 * @param {VectorXd} Xsig_aug a single sigma point
 */
void UKF::Predictionmodel(VectorXd Xsig_aug, double delta_t) {
	//extract values for better readability
	double p_x = Xsig_aug(0);
	double p_y = Xsig_aug(1);
	double v = Xsig_aug(2);
	double yaw = Xsig_aug(3);
	double yawd = Xsig_aug(4);
	double nu_a = Xsig_aug(5);
	double nu_yawdd = Xsig_aug(6);

	//predicted state values
	double px_p, py_p;

	//avoid division by zero
	if (fabs(yawd) > 0.001) {
		px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
		py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
	}
	else {
		px_p = p_x + v * delta_t*cos(yaw);
		py_p = p_y + v * delta_t*sin(yaw);
	}

	double v_p = v;
	double yaw_p = yaw + yawd * delta_t;
	double yawd_p = yawd;

	//add noise
	px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
	py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
	v_p = v_p + nu_a * delta_t;

	yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
	yawd_p = yawd_p + nu_yawdd * delta_t;

	//write predicted sigma point into right column
	Xsig_pred_(0, i) = px_p;
	Xsig_pred_(1, i) = py_p;
	Xsig_pred_(2, i) = v_p;
	Xsig_pred_(3, i) = yaw_p;
	Xsig_pred_(4, i) = yawd_p;
}
 
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	/**
	Estimate the object's location. Modify the state vector, x_. Predict sigma points, 
	the state, and the state covariance matrix.
	*/

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//create augmantaion matrix
	augmantaion(Xsig_aug);

	//predict sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		Predictionmodel(Xsig_aug.col(i), delta_t);
	}

	//predicted state mean
	x_.fill(0.0);

	x_ = x_ + firstweight_ * Xsig_pred_.col(0);
	for (int i = 1; i < 2 * n_aug_ + 1; i++) 
	{ 
		//iterate over sigma points
		x_ = x_ + remainingweights_ * Xsig_pred_.col(i);
	}

	//predicted state covariance matrix
	P_.fill(0.0);
	VectorXd x_diff = Xsig_pred_.col(0) - x_;
	P_ = P_ + firstweight_ * x_diff * x_diff.transpose();
	for (int i = i; i < 2 * n_aug_ + 1; i++)
	{  //iterate over sigma points

		// state difference
		x_diff = Xsig_pred_.col(i) - x_;
		P_ = P_ + remainingweights_ * x_diff * x_diff.transpose();
	}

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	Use lidar data to update the belief about the object's position. 
	Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/
	int n_z_ = 2;
	VectorXd x_diff, z_diff;
	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_);
	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z_, n_z_);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

	  // extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);

		// measurement model
		Zsig(0, i) = p_x;                        //x
		Zsig(1, i) = p_y;						//y
	}

	z_pred.fill(0.0);
	z_pred = z_pred + firstweight_ * Zsig.col(0);
	for (int i = 1; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + remainingweights_ * Zsig.col(i);
	}

	//innovation covariance matrix S
	MatrixXd S = MatrixXd(n_z_, n_z_);
	S.fill(0.0);
	z_diff = Zsig.col(0) - z_pred;
	S = S +firstweight_ * z_diff * z_diff.transpose();

	for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		z_diff = Zsig.col(i) - z_pred;
		S = S + remainingweights_ * z_diff * z_diff.transpose();
	}

	R << std_laspx_ * std_laspx_, 0,
		0, std_laspy_*std_laspy_;
	S = S + R;

	//calculate cross correlation matrix
	Tc.fill(0.0);

	//residual
	z_diff = Zsig.col(0) - z_pred;

	// state difference
	x_diff = Xsig_pred_.col(0) - x_;

	//angle normalization
	Tc = Tc + firstweight_ * x_diff * z_diff.transpose();
	for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		//residual
		z_diff = Zsig.col(i) - z_pred;

		// state difference
		x_diff = Xsig_pred_.col(i) - x_;

		//angle normalization
		Tc = Tc + remainingweights_ * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	z_diff = meas_package.raw_measurements_ - z_pred;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;

	P_ = P_ - K * S*K.transpose();

	NIS_lidar.push_back(z_diff.transpose() * S.inverse() * z_diff);

#ifdef DEBUG
	ChiSquare_FIVEpec(meas_package, CHI_lidar);
#endif // DEBUG
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/
	int n_z_ = 3;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_);
	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z_, n_z_);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

	  // extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0, i) = sqrt(p_x*p_x + p_y * p_y);                        //r
		Zsig(1, i) = atan2(p_y, p_x);                                 //phi
		Zsig(2, i) = (p_x*v1 + p_y * v2) / sqrt(p_x*p_x + p_y * p_y);   //r_dot
	}


	z_pred.fill(0.0);
	z_pred = z_pred + firstweight_ * Zsig.col(i);
	for (int i = 1; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + remainingweights_ * Zsig.col(i);
	}

	//innovation covariance matrix S
	MatrixXd S = MatrixXd(n_z_, n_z_);
	S.fill(0.0);

	VectorXd z_diff = Zsig.col(0) - z_pred;
	//angle normalization
	while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
	S = S + firstweight_ * z_diff * z_diff.transpose();
	for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

		S = S + remainingweights_ * z_diff * z_diff.transpose();
	}

	R << std_radr_ * std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;
	S = S + R;

	//calculate cross correlation matrix
	Tc.fill(0.0);

	//residual
	VectorXd z_diff = Zsig.col(i) - z_pred;
	//angle normalization
	while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

	// state difference
	VectorXd x_diff = Xsig_pred_.col(i) - x_;

	//angle normalization
	while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
	while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
	
	Tc = Tc + firstweight_ * x_diff * z_diff.transpose();

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

	  //residual
		z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		x_diff = Xsig_pred_.col(i) - x_;

		//angle normalization
		while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + remainingweights_ * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

	//angle normalization
	while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S*K.transpose();

	NIS_radar.push_back(z_diff.transpose() * S.inverse() * z_diff);
#ifdef DEBUG
	ChiSquare_FIVEpec(meas_package, CHI_radar);
#endif
}

#ifdef DEBUG
void UKF::ChiSquare_FIVEpec(MeasurementPackage meas_package, double CHI_value)
{
	double percentage = 0;

	if (meas_package.sensor_type_ == MeasurementPackage::LASER)
	{

		for (int index = 0; index < NIS_lidar.size(); index++)
		{
			if (NIS_lidar[index] >= CHI_lidar)
			{
				percentage++;
			}

		}
		std::cout << "laser performance = " << (100 * percentage) / NIS_lidar.size() << std::endl;
	}
	else
	{
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
			for (int index = 0; index < NIS_radar.size(); index++)
			{
				if (NIS_radar[index] >= CHI_radar)
				{
					percentage++;
				}

			}
			std::cout << "radar performance = " << (100 * percentage) / NIS_radar.size() << std::endl;
		}
	}


}
#endif
