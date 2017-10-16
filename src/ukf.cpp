#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define ALMOST_ZERO 0.0001

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // State dimension includes px, py, v, psi, psi_dot
  n_x_ = x_.size();
  
  // Augmented state dimension includes px, py, v, psi, psi_dot, nu_a, nu_psi_dot_dot
  n_aug_ = n_x_ + 2;

  // Define sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Define the sigma points amount
  n_sigma_ = 2*n_aug_ + 1;

  // Create and define vector for weights
  weights_ = VectorXd(n_sigma_);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int index = 1; index < n_sigma_; index++) {
    weights_(index) = 1/(2*(lambda_ + n_aug_));
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "UKF: " << endl;
    x_ = VectorXd(n_x_);

    //state covariance matrix P
    P_ = MatrixXd::Identity(n_x_, n_x_);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to CTRV model and initialize state.
      */
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      double px = rho*sin(phi);
      double py = rho*cos(phi);
      double vx = rho_dot*sin(phi);
      double vy = rho_dot*cos(phi);
      double vx2_vy2 = vx*vx + vy*vy;
      double vx2_vy2_sqrt = sqrt(vx2_vy2);

      tools.KeepNoneZero(vx2_vy2_sqrt, ALMOST_ZERO);

      x_ << px, py, vx2_vy2/vx2_vy2_sqrt, phi, 0;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      //double phi = atan2(py, px);
      //x_ << px, py, 0, phi, 0;
      x_ << px, py, 0, 0, 0;

    }

    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  float dt = (meas_package.timestamp_ - time_us_)/1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  
  if (dt > ALMOST_ZERO) {
    Prediction(dt);
  }
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar update
    if (use_radar_) {
      UpdateRadar(meas_package);
    }
  } else {
    if (use_laser_) {
      UpdateLidar(meas_package);
    }
  }
}  
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /** 
   * First, create augmented sigma points
  */

  // Create augmented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);

  // Create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);

  // Create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, n_sigma_); 

  // Define the augmented mean state
  for (int i = 0; i < n_x_; i++) {
      x_aug_(i) = x_(i);
  }
  x_aug_(n_x_) = 0;
  x_aug_(n_x_ + 1) = 0;

  // Create and define augmented covariance matrix
  MatrixXd Q_ = MatrixXd(2, 2);
  Q_ << std_a_*std_a_, 0,
       0, std_yawdd_*std_yawdd_;

  // Define augmented state covariance matrix
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_.bottomRightCorner(2, 2) = Q_;
  P_aug_.topRightCorner(n_x_, 2) = MatrixXd::Zero(n_x_, 2);
  P_aug_.bottomLeftCorner(2, n_x_) = MatrixXd::Zero(2, n_x_);

  // Create square root matrix
  MatrixXd A_ = P_aug_.llt().matrixL();

  // Define augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;
  float lambda_naug = sqrt(lambda_ + n_aug_);
  VectorXd A_l = VectorXd(n_aug_);
  int col_num;
  for (int i = 0; i < n_aug_; i++)
  {
    A_l = lambda_naug*A_.col(i);
    col_num = i + 1;
    Xsig_aug_.col(col_num) = x_aug_ + A_l;
    Xsig_aug_.col(col_num + n_aug_) = x_aug_ - A_l;
  }

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);
  float delta_t_2 = delta_t*delta_t;
  for (int col_num = 0; col_num < n_sigma_; col_num++) {
      VectorXd X_col = Xsig_aug_.col(col_num);
      float v = X_col(2);
      float phi = X_col(3);
      float phi_dot = X_col(4);
      float nu_a = X_col(5);
      float nu_phi_dot_dot = X_col(6);
      VectorXd main_part = VectorXd(n_x_);
      VectorXd nu_influence = VectorXd(n_x_);
      if (fabs(phi_dot) < ALMOST_ZERO) {
          main_part << v*cos(phi)*delta_t,
                       v*sin(phi)*delta_t,
                       0,
                       phi_dot*delta_t,
                       0;
          
      } else {
          main_part << v*(sin(phi+phi_dot*delta_t)-sin(phi))/phi_dot,
                       v*(-cos(phi+phi_dot*delta_t)+cos(phi))/phi_dot,
                       0,
                       phi_dot*delta_t,
                       0;
          
      }
      nu_influence << delta_t_2*cos(phi)*nu_a/2,
                      delta_t_2*sin(phi)*nu_a/2,
                      delta_t*nu_a,
                      delta_t_2*nu_phi_dot_dot/2,
                      delta_t*nu_phi_dot_dot;

      Xsig_pred_.col(col_num) = X_col.head(n_x_) + main_part + nu_influence;
  }

  // Predict state mean
  x_ = Xsig_pred_*weights_;
  
  // Predict state covariance matrix
  P_ = MatrixXd::Zero(n_x_, n_x_);
  for (int index = 0; index < n_sigma_; index++) {
      MatrixXd diff = Xsig_pred_.col(index) - x_;
      P_ = P_ + weights_(index)*diff*diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  //transform sigma points into measurement space
  Zsig = Xsig_pred_.block(0, 0, 2, n_sigma_); 
 
  //calculate mean predicted measurement
  z_pred = Zsig*weights_;

  // Create temperary vector to calculate measurement covariance matrix S
  VectorXd Zsig_diff = VectorXd(n_z);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // Create temperary vector to calculate cross correlation matrix
  VectorXd Xsig_diff = VectorXd(n_x_);
 
  for (int index = 0; index < n_sigma_; index++) {
      Zsig_diff = Zsig.col(index) - z_pred;
      S = S + weights_(index)*Zsig_diff*Zsig_diff.transpose();
      Xsig_diff = Xsig_pred_.col(index) - x_;
      Tc = Tc + weights_(index)*Xsig_diff*Zsig_diff.transpose();  //matrix 5x2
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  
  S = S + R;
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_;

  //calculate Kalman gain K;
  MatrixXd S_inverse = S.inverse();
  MatrixXd K = Tc*S_inverse;
  
  //update state mean and covariance matrix
  VectorXd z_dif = z - z_pred;
  x_ = x_ + K*z_dif;
  P_ = P_ - K*S*K.transpose();

  // Calculate the NIS
  NIS_laser_ = z_dif.transpose()*S_inverse*z_dif; 

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
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  //transform sigma points into measurement space
  double px, py, v, phi, phi_dot, px2_py2_sqrt;
  VectorXd Xsig_pred_col = VectorXd(n_x_);
  VectorXd Zsig_col = VectorXd(n_z);
  for (int index = 0; index < n_sigma_; index++) {
      Xsig_pred_col = Xsig_pred_.col(index);
      px = Xsig_pred_col(0);
      py = Xsig_pred_col(1);
      v = Xsig_pred_col(2);
      phi = Xsig_pred_col(3);
      phi_dot = Xsig_pred_col(4);
      px2_py2_sqrt = sqrt(px*px + py*py);
      Zsig_col << px2_py2_sqrt,
                  atan2(py, px),
                  (px*cos(phi)*v + py*sin(phi)*v)/px2_py2_sqrt;
      Zsig.col(index) = Zsig_col;
  }
  
  //calculate mean predicted measurement
  z_pred = Zsig*weights_;
  
  //Create temperary vector to calculate measurement covariance matrix S
  VectorXd Zsig_diff = VectorXd(n_z);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // Create temperary vector to calculate cross correlation matrix
  VectorXd Xsig_diff = VectorXd(n_x_);
 
  for (int index = 0; index < n_sigma_; index++) {
      Zsig_diff = Zsig.col(index) - z_pred;
      S = S + weights_(index)*Zsig_diff*Zsig_diff.transpose();
      Xsig_diff = Xsig_pred_.col(index) - x_;
      Tc = Tc + weights_(index)*Xsig_diff*Zsig_diff.transpose();  //matrix 5x3
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0 , 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  S = S + R;
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_;

  //calculate Kalman gain K;
  MatrixXd S_inverse = S.inverse();
  MatrixXd K = Tc*S.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_dif = z - z_pred;
  tools.KeepDiffInTwoPi(z_dif(1), M_PI);
  x_ = x_ + K*z_dif;
  P_ = P_ - K*S*K.transpose();
  
  // Calculate the NIS
  NIS_radar_ = z_dif.transpose()*S_inverse*z_dif; 

}
