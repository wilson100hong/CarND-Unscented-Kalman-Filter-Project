#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse.fill(0);
  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    cout << "CalculateRMSE(): Invalid estimations or ground_truth" << endl;
    return rmse;
  }
  for (size_t i=0;i<estimations.size();++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}
