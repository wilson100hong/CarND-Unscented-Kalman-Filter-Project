#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE here.
   */
  int n = estimations.size();
  VectorXd rmse(n);
  for (int i=0;i<n;++i) {
    rsme[i] = 0;
  }
  if (n == 0 || n != ground_truth.size()) {
    cout << "CalculateRMSE(): Invalid estimations or ground_truth" << endl;
    return rmse;
  }
  
  for (size_t i=0;i<n;++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  rmse = rmse / n;
  rmse = rmse.array().sqrt();
  return rmse;
}
