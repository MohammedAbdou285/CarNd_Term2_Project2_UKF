#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
  {
      VectorXd rmse(4);
      rmse << 0, 0, 0, 0;

      // check if there is a mismatch sizes between the estimations and the groundtruth
      //		 if the estimations size is not equal to zero		
      if (estimations.size() != ground_truth.size() || estimations.size() == 0)
      {
        cout << "Mismatch Sizes" << endl;
      }

      else
      {
        // accumulate squared residuals
        for (unsigned int i = 0; i < estimations.size(); ++i)
        {
          // Define the residual Vector
          VectorXd residual = estimations[i] - ground_truth[i];
          // Coefficient-wise multiplication
          residual = residual.array() * residual.array();
          rmse += residual;
        }

        // Calculate mean
        rmse = rmse / estimations.size();
        // Calculate the square root of the mean to have the root mean squared
        rmse = rmse.array().sqrt();
      }

      // return zeros in case of mismatch sizes or 
      // return the calculated RMSE between the estimations and the groundtruth 
      return rmse;
  }