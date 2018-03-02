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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
  
  // ************************TODO*************************************
  ///* set initialization flag to be false
  is_initialized_ = false;

  ///* set time to be zero initially
  time_us_ = 0;
  
  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_+ 1);

  ///* NIS_radar
  NIS_radar_ = 0.0;

  ///* NIS_radar
  NIS_lidar_ = 0.0;




  //*******************************************************************
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

  /**********************************************************************
   *    Initialization
   *********************************************************************/

  if (!is_initialized_)
  {
    x_ << 1, 1, 1, 1, 0.1;
    P_ = MatrixXd::Identity(5,5);
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    {
        // Convert radar measurement from polar to cartesian coordinates
        float rho = meas_package.raw_measurements_(0);
        float phi = meas_package.raw_measurements_(1);
        float rho_dot = meas_package.raw_measurements_(2);

        float x = rho * cos(phi);
        float y = rho * sin(phi);

        x_(0) = x;
        x_(1) = y;
    }

    else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    {
        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    // It is only initialization, so no need to predict or update
    return;
  }

  /**********************************************************************
   *    Prediction
   *********************************************************************/

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  Prediction(dt);
  time_us_ = meas_package.timestamp_;

  
  /**********************************************************************
   *    Update
   *********************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
      UpdateRadar(meas_package);
  }

  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
      UpdateLidar(meas_package);
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

  /**********************************************************************
   *    1- Create or Generate Sigma points
   *********************************************************************/
  
  // Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // Create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // Create sigma points matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Create augmenetd covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  // Create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Create augmenetd sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++)
  {
      Xsig_aug.col(i+1)          = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
      Xsig_aug.col(i+1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }


  /**********************************************************************
   *    2- Predict Sigma Points
   *********************************************************************/

  for (int i = 0; i < 2*n_aug_+1; i++)
  {
      // Extract values for readability
      double p_x      = Xsig_aug(0,i); 
      double p_y      = Xsig_aug(1,i);
      double v        = Xsig_aug(2,i);
      double yaw      = Xsig_aug(3,i);
      double yawd     = Xsig_aug(4,i);
      double nu_a     = Xsig_aug(5,i);
      double nu_yawdd = Xsig_aug(6,i);

      // predicted state values
      double px_p, py_p;

      // avoid divison ny zero 
      if (fabs(yawd) > 0.001)
      {
          // General equations
          px_p = p_x + (v/yawd) * ( sin(yaw + yawd * delta_t) - sin(yaw));
          py_p = p_y + (v/yawd) * (-cos(yaw + yawd * delta_t) + cos(yaw));
      }
      else
      {
          // Special case
          px_p = p_x + v * delta_t * cos(yaw);
          py_p = p_y + v * delta_t * sin(yaw);
      }

      double v_p = v;
      double yaw_p = yaw + yawd * delta_t;
      double yawd_p = yawd;

      // add noise
      double dt2 = delta_t * delta_t;
      px_p = px_p + 0.5 * nu_a * dt2 * cos(yaw);
      py_p = py_p + 0.5 * nu_a * dt2 * sin(yaw);
      v_p = v_p + nu_a * delta_t;

      yaw_p  = yaw_p + 0.5 * nu_yawdd * dt2;
      yawd_p = yawd_p + nu_yawdd * delta_t;

      // write predicted sigma point into right column
      Xsig_pred_(0,i) = px_p; 
      Xsig_pred_(1,i) = py_p;
      Xsig_pred_(2,i) = v_p;
      Xsig_pred_(3,i) = yaw_p;
      Xsig_pred_(4,i) = yawd_p;
  }

  /**********************************************************************
   *    3- Calculate New Mean Vector and Covriance Matrix
   *********************************************************************/

  // set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  // calculate other weights 
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
      double weight = 0.5 / (n_aug_ + lambda_);
      weights_(i) = weight;
  }

  // predicted state mean
  x_.fill(0.0);
  //iterate over the all sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
      x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      // angle normalization
      while (x_diff(3) > M_PI)
          x_diff(3) -= 2.0 * M_PI;
      while (x_diff(3) < -M_PI)
          x_diff(3) += 2.0 * M_PI;

      P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
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
}
