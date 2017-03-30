#include <iostream>
#include "ukf.h"
#include "tools.h"

using namespace std;
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
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.9;

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
  n_x_ = 5;

  n_aug_ = 7;

  n_z_laser_ = 2;

  n_z_radar_ = 3;


  lambda_ = 3 - n_aug_ ;

  is_initialized_ = false;

  previous_timestamp_ = 0;

  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_+1);

  weights = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for(int i = 1; i<2*n_aug_+1;i++){
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
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
  if(!is_initialized_){
    std::cout<<"Kalman Filter Initialization"<<std::endl;
    x_ = VectorXd(n_x_);
    x_ << 0,0,0,0,0;

    P_ = MatrixXd(n_x_,n_x_);

    P_ <<  0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;


    time_us_ = meas_package.timestamp_;

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){


    float rho = meas_package.raw_measurements_[0];
    float phi = meas_package.raw_measurements_[1];
    float rhod = meas_package.raw_measurements_[2];

    float x_cart = rho * cos(phi);
    float y_cart = rho * sin(phi);
    float v = rhod;

    if(x_cart == 0 or y_cart == 0){
      return;
    }

    x_ << x_cart,y_cart,v,phi,0;

  }else if(meas_package.sensor_type_==MeasurementPackage::LASER){

    if(meas_package.raw_measurements_[0]==0 or meas_package.raw_measurements_[1]==0){
      return;
    }

    x_<<meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],0,0,0;
  }

  //done initializing, no need to predict or update
  is_initialized_ = true;
  std::cout<<"Initialized"<<std::endl;
  return;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds


  //Prediction
  while(dt>0.1){
    Prediction(0.05);
    dt-=0.05;
  }


  Prediction(dt);
  time_us_ = meas_package.timestamp_;



  //Measurement Update(RADAR)
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    UpdateRadar(Xsig_pred_,meas_package);



  //Measurement Update(LIDAR)
  }else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

    UpdateLidar(Xsig_pred_,meas_package);

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
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_,2*n_aug_+1);
  GenerateAugmentedSigmaPoints(x_,P_,Xsig_aug_);
  SigmaPointPrediction(Xsig_aug_,delta_t,Xsig_pred_);
  PredictMeanAndCovariance(Xsig_pred_,x_,P_);
}

 //Generating Augmented Sigma Points
 void UKF::GenerateAugmentedSigmaPoints(const VectorXd& x_, const MatrixXd& P_, MatrixXd& Xsig_aug_){
  lambda_ = 3 - n_aug_;
  VectorXd x_aug_ = VectorXd(n_aug_);
  MatrixXd P_aug_ = MatrixXd(n_aug_,n_aug_);

  x_aug_.head(5) = x_;
  x_aug_(5)= 0;
  x_aug_(6)=0;
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;
  MatrixXd L_ = P_aug_.llt().matrixL();
  Xsig_aug_.col(0) = x_aug_;

  for(int i = 0;i<n_aug_;i++){
    Xsig_aug_.col(i+1) = x_aug_ + sqrt(lambda_+n_aug_)*L_.col(i);
    Xsig_aug_.col(n_aug_+i+1) = x_aug_ - sqrt(lambda_+n_aug_)*L_.col(i);
  }

 }

 //Sigma Points Prediction
 void UKF::SigmaPointPrediction(const MatrixXd& Xsig_aug_, double delta_t, MatrixXd& Xsig_pred_){

  for(int i = 0 ; i < 2*n_aug_+1 ; i++){
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a =Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);


    double px_p,py_p;
    if(fabs(yawd)>0.001){
      px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );

    }else{
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);

    }
    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //transformation
    px_p = px_p + 0.5*nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;
    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;

    }

 }

 //Calculate Predicted Mean and Covariance
 void UKF::PredictMeanAndCovariance(const MatrixXd& Xsig_pred_, VectorXd& x_, MatrixXd& P_){

    x_.fill(0.0);
    for(int i=0; i<2*n_aug_+1;i++){
      x_ = x_ + weights(i)*Xsig_pred_.col(i);
    }
    P_.fill(0.0);
    for(int i = 0 ; i<2*n_aug_+1;i++){
      VectorXd x_diff = Xsig_pred_.col(i)-x_;
      while(x_diff(3)>M_PI) x_diff(3) -=2.*M_PI;
      while(x_diff(3)<-M_PI) x_diff(3) +=2.*M_PI;

    P_ = P_ + weights(i) * x_diff*x_diff.transpose();
    }

 }

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MatrixXd& Xsig_pred_ ,MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  VectorXd z_pred = VectorXd(n_z_laser_);
  MatrixXd S_laser = MatrixXd(n_z_laser_,n_z_laser_);

  MatrixXd Zsig = MatrixXd(n_z_laser_,2*n_aug_+1);
  PredictLidarMeasurement(Xsig_pred_, z_pred, S_laser, Zsig);
  VectorXd z_ = meas_package.raw_measurements_;
  UpdateLidarState(z_,Xsig_pred_,z_pred,Zsig,S_laser,x_,P_);
}

void UKF::PredictLidarMeasurement(const MatrixXd& Xsig_pred_, VectorXd& z_pred, MatrixXd& S_laser, MatrixXd& Zsig){

  for(int i = 0; i<2*n_aug_+1;i++){
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);


    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }


  z_pred.fill(0.0);
  for(int i = 0; i<2*n_aug_+1;i++){
    z_pred = z_pred + weights(i)*Zsig.col(i);
  }


  S_laser.fill(0.0);
  for(int i = 0; i<2*n_aug_+1;i++){
    VectorXd z_diff = Zsig.col(i)-z_pred;
    S_laser = S_laser + weights(i)*z_diff*z_diff.transpose();
  }

  MatrixXd R_laser = MatrixXd(n_z_laser_,n_z_laser_);
  R_laser << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;

  S_laser = S_laser + R_laser;
}

void UKF::UpdateLidarState(const VectorXd& z_,const MatrixXd& Xsig_pred_, const VectorXd& z_pred, const MatrixXd& Zsig,
                      const MatrixXd& S_laser, VectorXd& x_ , MatrixXd& P_){
  MatrixXd Tc_laser = MatrixXd(n_x_,n_z_laser_);
  Tc_laser.fill(0.0);
  for(int i = 0; i<2*n_aug_+1;i++){
    VectorXd z_diff = Zsig.col(i)-z_pred;
    VectorXd x_diff = Xsig_pred_.col(i)-x_;

    Tc_laser = Tc_laser + weights(i)*x_diff*z_diff.transpose();
  }

  MatrixXd K_laser = Tc_laser * S_laser.inverse();
  VectorXd z_diff = z_-z_pred;
  x_ = x_ + K_laser*z_diff;
  P_ = P_ - K_laser*S_laser*K_laser.transpose();

  NIS_laser_ = z_diff.transpose()*S_laser.inverse()*z_diff;


  //cout<<x_<<endl;
  //cout<<P_<<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MatrixXd& Xsig_pred_,MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  VectorXd z_pred = VectorXd(n_z_radar_);
  MatrixXd S_radar = MatrixXd(n_z_radar_,n_z_radar_);
  MatrixXd Zsig = MatrixXd(n_z_radar_,2*n_aug_+1);
  PredictRadarMeasurement(Xsig_pred_, z_pred, S_radar, Zsig);
  VectorXd z_ = meas_package.raw_measurements_;
  UpdateRadarState(z_,Xsig_pred_,z_pred,Zsig,S_radar,x_,P_);
}

void UKF::PredictRadarMeasurement(const MatrixXd& Xsig_pred_, VectorXd& z_pred, MatrixXd& S_radar, MatrixXd& Zsig){

  for (int i = 0; i<2*n_aug_+1;i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double vx = cos(yaw) * v;
    double vy = sin(yaw) * v;

    Zsig(0,i) = sqrt(p_x*p_x+p_y*p_y);

    Zsig(1,i) = atan2(p_y,p_x);
    if(Zsig(0,i) > 0.001){
      Zsig(2,i) = (p_x*vx + p_y*vy)/sqrt(p_x*p_x + p_y*p_y);

    }else{

      Zsig(2,i) = 0.0;
    }

  }

  z_pred.fill(0.0);
  for(int i = 0; i<2*n_aug_+1;i++){
    z_pred = z_pred + weights(i)*Zsig.col(i);
  }

  S_radar.fill(0.0);
  for(int i = 0; i<2*n_aug_+1;i++){
    VectorXd z_diff = Zsig.col(i)- z_pred;
    while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S_radar = S_radar + weights(i)*z_diff*z_diff.transpose();
  }

  MatrixXd R_radar = MatrixXd(n_z_radar_,n_z_radar_);
  R_radar << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;
  S_radar = S_radar + R_radar;

}

void UKF::UpdateRadarState(const VectorXd& z_,const MatrixXd& Xsig_pred_, const VectorXd& z_pred, const MatrixXd& Zsig,
                      const MatrixXd& S_radar, VectorXd& x_ , MatrixXd& P_){
  MatrixXd Tc_radar = MatrixXd(n_x_,n_z_radar_);
  Tc_radar.fill(0.0);
  for(int i = 0 ; i<2*n_aug_+1; i++){
    //redidual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    //state difference
    VectorXd x_diff = Xsig_pred_.col(i)-x_;
    while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
    while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc_radar = Tc_radar + weights(i)*x_diff*z_diff.transpose();
  }
    MatrixXd K_radar = Tc_radar * S_radar.inverse();
    VectorXd z_diff =  z_-z_pred;
    while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    x_ = x_ + K_radar*z_diff;
    P_ = P_ - K_radar*S_radar*K_radar.transpose();

    NIS_radar_ = z_diff.transpose()*S_radar.inverse()*z_diff;


    //cout<<x_<<endl;
    //cout<<P_<<endl;
}