#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  *  A helper method to keep a variable none zero value
  */
  void KeepNoneZero(double &x, const float almost_zero);

  /**
  * Keep the angle difference in 2*PI
  * @param phi_diff The angle difference
  */
  void KeepDiffInTwoPi(double &phi_diff, const float pi);

};

#endif /* TOOLS_H_ */
