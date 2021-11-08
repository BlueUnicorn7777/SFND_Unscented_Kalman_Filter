#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

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
    std_a_ = 30;

    // Process noise standard deviation yaw acceleration in rad/s^2
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
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
    is_initialized_ = false;
    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3-n_aug_;
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    weights_ = VectorXd(2*n_aug_+1);
    time_us_ = 0;

    std_a_ = 2;
    std_yawdd_ = 1;


    pred_.x_aug =VectorXd(n_aug_);
    pred_.P_aug=MatrixXd(n_aug_, n_aug_);
    pred_.Xsig_aug=MatrixXd(n_aug_, 2 * n_aug_ + 1);
    pred_.x_aug.fill(0);
    pred_.P_aug.fill(0);
    //augumented process noise
    pred_.P_aug(5,5) = std_a_*std_a_;
    pred_.P_aug(6,6) = std_yawdd_*std_yawdd_;

    pred_.A =lambda_+n_aug_;


    weights_(0) = lambda_/(lambda_+n_aug_);
    double weight = 0.5/(n_aug_+lambda_);
    for (int i=1; i<2*n_aug_+1; ++i) {  // 2n+1 weights
        weights_(i) = weight;
    }
    r.n_z=3;
    //radar measurement noise covariance
    r.R = MatrixXd(r.n_z,r.n_z);
    r.R <<  std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0,std_radrd_*std_radrd_;
    r.Zsig = MatrixXd(r.n_z, 2 * n_aug_ + 1);
    r.z_pred = VectorXd(r.n_z);
    r.S = MatrixXd(r.n_z,r.n_z);
    r.Tc = MatrixXd(n_x_, r.n_z);
    l.K = MatrixXd(n_x_, r.n_z);

    l.n_z=2;
     //lidar measurement noise covariance
    l.R = MatrixXd(l.n_z,l.n_z);
    l.R.fill(0);
    l.R(0, 0) = std_laspx_*std_laspx_;
    l.R(1, 1) = std_laspy_*std_laspy_;

    l.Zsig = MatrixXd(l.n_z, 2 * n_aug_ + 1);
    l.z_pred = VectorXd(l.n_z);
    l.S = MatrixXd(l.n_z,l.n_z);
    l.Tc = MatrixXd(n_x_, l.n_z);
    l.K = MatrixXd(n_x_, l.n_z);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
    if (!is_initialized_) {
        if (meas_package.sensor_type_ == MeasurementPackage::LASER ) {
            x_ << meas_package.raw_measurements_[0],
                    meas_package.raw_measurements_[1],
                    0, 0 , 0;
            P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
                    0, std_laspy_*std_laspy_, 0, 0, 0,
                    0, 0, 1, 0, 0,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1;

        }
        else
            if (meas_package.sensor_type_ == MeasurementPackage::RADAR ) {
                double rho = meas_package.raw_measurements_(0);
                double phi = meas_package.raw_measurements_(1);
                double rhodot = meas_package.raw_measurements_(2);
                double x = rho * cos(phi);
                double y = rho * sin(phi);
                double vx = rhodot * cos(phi);
                double vy = rhodot * sin(phi);
                double v = sqrt(vx * vx + vy * vy);
                x_ << x, y, v, rho, rhodot;
                P_ << std_radr_*std_radr_, 0, 0, 0, 0,
                        0, std_radr_*std_radr_, 0, 0, 0,
                        0, 0, std_radrd_*std_radrd_, 0, 0,
                        0, 0, 0, std_radphi_, 0,
                        0, 0, 0, 0, std_radphi_;
            }

        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }

    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    //Predicition
     Prediction(dt);
    // Measurement update

    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
        //Prediction(dt);
        UpdateLidar(meas_package);}
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
        //Prediction(dt);
        UpdateRadar(meas_package);}

//    std::cout<<meas_package.raw_measurements_.x()<<"\t"<<meas_package.raw_measurements_.y()
//            ;

}

void UKF::Prediction(double delta_t) {
    /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
    //create augumented sigma points

    pred_.x_aug.head(5) = x_;
    pred_. P_aug.topLeftCorner(5,5) = P_;
    MatrixXd L = pred_.P_aug.llt().matrixL();
    pred_.Xsig_aug.col(0)  = pred_.x_aug;

    for (int i = 0; i< n_aug_; ++i) {
        pred_.Xsig_aug.col(i+1)       = pred_.x_aug + sqrt(pred_.A) * L.col(i);
        pred_.Xsig_aug.col(i+1+n_aug_) = pred_.x_aug - sqrt(pred_.A) * L.col(i);
    }

    // predict sigma points
    for (int i = 0; i< 2*n_aug_+1; ++i) {
        // extract values for better readability
        double p_x = pred_.Xsig_aug(0,i);
        double p_y = pred_.Xsig_aug(1,i);
        double v = pred_.Xsig_aug(2,i);
        double yaw = pred_.Xsig_aug(3,i);
        double yawd = pred_.Xsig_aug(4,i);
        double nu_a = pred_.Xsig_aug(5,i);
        double nu_yawdd = pred_.Xsig_aug(6,i);
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

    // predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }

    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
    }

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
    // transform sigma points into measurement space
    for (int i=0; i<2 * n_aug_ + 1; i++){
        l.Zsig(0, i) = Xsig_pred_(0, i);
        l.Zsig(1, i) = Xsig_pred_(1, i);
    }

    // calculate mean predicted measurement
    l.z_pred.fill(0);
    for (int i=0; i<2 * n_aug_ + 1; i++){
        l.z_pred += weights_(i) * l.Zsig.col(i);
    }
    // predicated measurement covariance matrix S
    l.S.fill(0);
    for (int i=0; i<2 * n_aug_ + 1; i++){
        // residual
        VectorXd z_diff = l.Zsig.col(i) - l.z_pred;
        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        l.S = l.S + weights_(i) * z_diff * z_diff.transpose();
    }
    // add measurement noise covariance matrix
    l.S += l.R;

    // lidar incoming measurement
    VectorXd z = VectorXd(l.n_z);
    z = meas_package.raw_measurements_;

    // calculate cross correlation matrix
    l.Tc.fill(0);
    for (int i = 0; i<2 * n_aug_ + 1; i++){
        // residual
        VectorXd z_diff = l.Zsig.col(i) - l.z_pred;
        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        l.Tc = l.Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // calculate Kalman gain K;
    l.K = l.Tc*l.S.inverse();

    // residual
    VectorXd z_diff = z - l.z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    //calculate NIS
     NIS_laser_ = z_diff.transpose() * l.S.inverse() * z_diff;
     //std::cout<<NIS_radar_<<std::endl;

    // update state mean and covariance matrix
   // x_ = x_ + l.K*(z - l.z_pred);
    x_ = x_ + l.K*(z_diff);
    P_ = P_ - l.K*l.S*l.K.transpose();

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

    r.Zsig.fill(0.0);
    r.z_pred.fill(0.0);
    r.S.fill(0.0);

    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // extract values for better readability
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);

        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        // sigma points in measurement space
        r.Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
        r.Zsig(1,i) = atan2(p_y,p_x);                                // phi
        r.Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   // r_dot
    }

    // predict sigma points

    r.z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; ++i) {
        r.z_pred = r.z_pred + weights_(i) * r.Zsig.col(i);
    }

    // predicated measurement covariance matrix S
    r.S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // residual
        VectorXd z_diff = r.Zsig.col(i) - r.z_pred;

        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        r.S = r.S + weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    r.S = r.S + r.R;

    r.Tc.fill(0.0);
    // create matrix for cross correlation Tc

    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // residual
        VectorXd z_diff = r.Zsig.col(i) -r.z_pred;
        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        r.Tc = r.Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // Kalman gain K;
    r.K = r.Tc * r.S.inverse();
    VectorXd z = VectorXd(r.n_z);
    z << meas_package.raw_measurements_;

    // residual
    VectorXd z_diff = z - r.z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    //calculate NIS
    NIS_radar_ = z_diff.transpose() * r.S.inverse() * z_diff;
    //std::cout<<NIS_radar_<<std::endl;

    // update state mean and covariance matrix
    x_ = x_ + r.K * z_diff;
    P_ = P_ - r.K*r.S*r.K.transpose();

}
