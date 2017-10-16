#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse;
  int m = estimations[0].size();
  rmse.setZero(m);

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
  	|| estimations.size() == 0){
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }
 
  //accumulate squared residuals
  for(unsigned int i = 0; i < estimations.size(); ++i){

    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;

}

void Tools::KeepNoneZero(double &x, const float almost_zero) {
  if (fabs(x) < almost_zero) {
    if (x > 0) {
      x = almost_zero;
    } else {
      x = -almost_zero;
    }
  }
}

void Tools::KeepDiffInTwoPi(double &phi_diff, const float pi) {
  if (phi_diff > pi) { 
    phi_diff = phi_diff - 2*pi;
  } else if (phi_diff < -pi) {
    phi_diff = phi_diff + 2*pi;
  }
}
