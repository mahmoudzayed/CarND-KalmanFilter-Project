#include "kalman_filter.h"
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
kalman::kalman() {

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

	///* State dimension
	n_x_ = 5;

	P_ = MatrixXd(n_x_, n_x_);

	x_ = VectorXd(n_x_);

#ifdef DEBUG
	CHI_radar = 7.815;

	CHI_lidar = 5.991;
#endif
}

kalman::~kalman() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void kalman::init()
{

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * Estimate the object's location. Modify the state vector, x_. Predict sigma points,
 * the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void kalman::Prediction(double delta_t) {
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
 * Predicts sigma points, the state, and the state covariance matrix.
 * Estimate the object's location. Modify the state vector, x_. Predict sigma points,
 * the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */

void kalman::Update()
{
	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();
	MatrixXd S = (S + S')/2;
	
	MatrixXd z_diff = meas_package.raw_measurements_ - z_pred;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;

	P_ = P_ - K * S*K.transpose();

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * Estimate the object's location. Modify the state vector, x_. Predict sigma points,
 * the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */

void kalman::ellipsoidalGating()
{
	/*
	nction [z_ingate, meas_in_gate] = ellipsoidalGating(state_pred, z, measmodel, gating_size)
            %ELLIPSOIDALGATING performs ellipsoidal gating for a single
            %object 
            %INPUT:  z: measurements --- (measurement dimension) x (number
            %        of measurements) matrix 
            %        state_pred: a structure with two fields:
            %                   x: predicted object state mean --- (state
            %                   dimension) x 1 vector 
            %                   P: predicted object state covariance ---
            %                   (state dimension) x (state dimension)
            %                   matrix  
            %        measmodel: a structure specifies the measurement model
            %        parameters 
            %        gating_size: gating size --- scalar
            %OUTPUT: z_ingate: measurements in the gate --- (measurement
            %        dimension) x (number of measurements in the gate)
            %        matrix
            %        meas_in_gate: boolean vector indicating whether the
            %        corresponding measurement is in the gate or not ---
            %        (number of measurements) x 1
            
            % initialize the output z_ingate
            z_ingate = [];
            % initaize the output meas_in_gate
            meas_in_gate = false(size(z,2),1);
            %Measurement model Jacobian
            Hx = measmodel.H(state_pred.x);
            %Innovation covariance
            S = Hx*state_pred.P*Hx' + measmodel.R;
            %Make sure matrix S is positive definite
            S = (S+S')/2;
            %K = (state_pred.P*Hx')/S;
            %State update
           % state_upd.x = state_pred.x + K*(z - measmodel.h(state_pred.x));
            z_par = measmodel.h(state_pred.x)
            
            for n = 1:size(z,2)
                % calculate the difference between z_i and predicted state
                z_diff = z(:,n) - z_par;
                G_threshold = z_diff' * (inv(S)) * z_diff;
                if G_threshold <= gating_size
                    z_ingate = [z_ingate z(:,n)];
                    meas_in_gate(n) = true;
                end
            end
	*/
}

void kalman::momentMatching()
{
	/*
	%MOMENTMATCHING: approximate a Gaussian mixture density as a
			%single Gaussian using moment matching
			%INPUT: w: normalised weight of Gaussian components in
			%       logarithm domain --- (number of Gaussians) x 1 vector
			%       states: structure array of size (number of Gaussian
			%       components x 1), each structure has two fields
			%               x: means of Gaussian components --- (variable
			%               dimension) x 1 vector
			%               P: variances of Gaussian components ---
			%               (variable dimension) x (variable dimension) matrix
			%OUTPUT:state: a structure with two fields:
			%               x: approximated mean --- (variable dimension) x
			%               1 vector
			%               P: approximated covariance --- (variable
			%               dimension) x (variable dimension) matrix

			if length(w) == 1
				state = states;
				return;
			end
			w = exp(w);
			state.x = zeros(size(states(1).x));
			state.P = zeros(size(states(1).P));

			for index = 1:length(states)
				state.x = state.x + (w(index)*states(index).x);
			end

			for index = 1:length(states)
				statediff = state.x - states(index).x
				state.P = state.P + (w(index)*states(index).P) + (w(index)*statediff*statediff');
			end
	*/
}

void kalman::mixtureReduction()
{
	/*
	%MIXTUREREDUCTION: uses a greedy merging method to reduce the
            %number of Gaussian components for a Gaussian mixture density 
            %INPUT: w: normalised weight of Gaussian components in
            %       logarithmic scale --- (number of Gaussians) x 1 vector 
            %       states: structure array of size (number of Gaussian
            %       components x 1), each structure has two fields 
            %               x: means of Gaussian components --- (variable
            %               dimension) x (number of Gaussians) matrix 
            %               P: variances of Gaussian components ---
            %               (variable dimension) x (variable dimension) x
            %               (number of Gaussians) matrix  
            %       threshold: merging threshold --- scalar
            %INPUT: w_hat: normalised weight of Gaussian components in
            %       logarithmic scale after merging --- (number of
            %       Gaussians) x 1 vector  
            %       states_hat: structure array of size (number of Gaussian
            %       components after merging x 1), each structure has two
            %       fields  
            %               x: means of Gaussian components --- (variable
            %               dimension) x (number of Gaussians after
            %               merging) matrix  
            %               P: variances of Gaussian components ---
            %               (variable dimension) x (variable dimension) x
            %               (number of Gaussians after merging) matrix  
            
            if length(w) == 1
                w_hat = w;
                states_hat = states;
                return;
            end
            
            %Index set of components
            I = 1:length(states);
            el = 1;
            
            while ~isempty(I)
                Ij = [];
                %Find the component with the highest weight
                [~,j] = max(w);

                for i = I
                    temp = states(i).x-states(j).x;
                    val = diag(temp.'*(states(j).P\temp));
                    %Find other similar components in the sense of small
                    %Mahalanobis distance 
                    if val < threshold
                        Ij= [ Ij i ];
                    end
                end
                
                %Merge components by moment matching
                [temp,w_hat(el,1)] = normalizeLogWeights(w(Ij));
                states_hat(el,1) = GaussianDensity.momentMatching(temp, states(Ij));
                
                %Remove indices of merged components from index set
                I = setdiff(I,Ij);
                %Set a negative to make sure this component won't be
                %selected again 
                w(Ij,1) = log(eps);
                el = el+1;
            end
            
        end
        
    end

	*/
}
