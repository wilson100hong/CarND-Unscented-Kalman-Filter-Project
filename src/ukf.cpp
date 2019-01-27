#include "ukf.h"
#include "Eigen/Dense"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

double normalize_rad(double rad) {
  while (rad > M_PI) rad -= 2*M_PI;
  while (rad < -M_PI) rad += 2*M_PI;
  return rad;
}

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
  // TODO: tuning
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // TODO: tuning
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(1 + 2*n_aug_);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1;i<2*n_aug_+1;++i) {
    weights_(i) = 0.5/(lambda_+n_aug_);
  }

  x_ = VectorXd(n_x_);
  x_.fill(0.0);
  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  Xsig_pred_.fill(0);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  auto sensor_type = meas_package.sensor_type_;
  if (!is_initialized_) {
    // Initalize state x_ and state covariance P_
    if (sensor_type == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      // and initialize state.
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      //double rho_dot = measurement_pack.raw_measurements_[2];
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      // TODO: make sure x_ is correctly set
      x_ << px, py, 0, 0, 0;
      // TODO: tune P_
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0, 
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    } else {
      // Initialize state from lidar. 
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      // TODO: make sure x_ is correctly set
      x_ << px, py, 0, 0, 0;
      // TODO: tune P_
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0, 
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    }
      
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }
  /**
   * Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  // Ignore sensor measurement if not enabled
  if ((sensor_type == MeasurementPackage::RADAR && !use_radar_) ||
      (sensor_type == MeasurementPackage::LASER && !use_laser_)) {
    return;
  }

  // Predict
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);

  // Update
  if (sensor_type == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }

  // TODO: debug
  cout << "x_ = " << x_ << endl; 
  cout << "P_ = " << P_ << endl; 
}

void UKF::Prediction(double delta_t) {
  /**
   * Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
   // Generate augmented sigma points: Xsig_aug
   // Augmented mean state.
   VectorXd x_aug = VectorXd(n_aug_);  
   x_aug.head(n_x_) = x_;
   x_aug(5) = 0;
   x_aug(6) = 0;
   
   // Augmented state covariance.
   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
   P_aug.fill(0.0);
   P_aug.topLeftCorner(n_x_, n_x_) = P_;
   P_aug(5, 5) = std_a_ * std_a_;
   P_aug(6, 6) = std_yawdd_* std_yawdd_;

   // Diagonal matrix of P_aug
   MatrixXd A = P_aug.llt().matrixL();

   // Augmented sigma points.
   MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
   Xsig_aug(0) = x_aug;
   for (int i=0;i<n_aug_;++i) {
     Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
     Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
   }

   // Predict sigma points
   Xsig_pred_.fill(0.0);
   for (int i=0;i<2*n_aug_+1;++i) {
     // Extract values for better readability
     double p_x = Xsig_aug(0,i);
     double p_y = Xsig_aug(1,i);
     double v = Xsig_aug(2,i);
     double yaw = Xsig_aug(3,i);
     double yawd = Xsig_aug(4,i);
     double nu_a = Xsig_aug(5,i);
     double nu_yawdd = Xsig_aug(6,i);

     // Predicted state values
     double px_p, py_p;

     // Avoid division by zero
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

     // Add noise
     px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
     py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
     v_p = v_p + nu_a*delta_t;

     yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
     yawd_p = yawd_p + nu_yawdd*delta_t;

     // Write predicted sigma point into right column
     Xsig_pred_(0,i) = px_p;
     Xsig_pred_(1,i) = py_p;
     Xsig_pred_(2,i) = v_p;
     Xsig_pred_(3,i) = yaw_p;
     Xsig_pred_(4,i) = yawd_p;
   }

   // Predict state mean
   x_.fill(0);
   for (int i=0;i<2*n_aug_+1;++i) {
     x_ = x_ + weights_(i)*Xsig_pred_.col(i);
   }
   // Predict state covariance
   P_.fill(0);
   for (int i=0;i<2*n_aug_+1;++i) {
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     // Need to normalize angle
     x_diff(3) = normalize_rad(x_diff(3));
     P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
   }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z = 2;
  MatrixXd H(n_z, n_x_);
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;
  // Measurement residual
  VectorXd y = meas_package.raw_measurements_ - H * x_;
  MatrixXd R(2, 2);
  R << std_laspx_* std_laspx_, 0,
       0, std_laspy_ * std_laspy_;

  MatrixXd S = H * P_ * H.transpose() + R;
  MatrixXd K = P_* H.transpose() * S.inverse();

  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P = (I - K * H) * P_;

  // Compute NIS
  double NIS = y.transpose() * S.inverse() * y;
  cout << "Radar NIS: " << NIS << endk;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z = 3;
  // Sigma points in measuresument space.
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  
  // Transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    // Extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // Measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   // r_dot
  }

  // Predicted measurement mean.
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Measurement covariance.
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // Need to normalize angle
    z_diff(1) = normalize_rad(z_diff(1);
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
  S = S + R;

  // Cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = normalize_rad(z_diff(1));

    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = normalize_rad(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // Measurement residual
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  y(1) = normalize_rad(y(1);

  // Update state mean and covariance matrix
  x_ = x_ + K * y;
  P_ = P_ - K * S * K.transpose();

  // Compute NIS
  double NIS = y.transpose() * S.inverse() * y;
  cout << "Radar NIS: " << NIS << endk;
}
