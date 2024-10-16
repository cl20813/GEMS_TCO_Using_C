#ifndef MATERN_COV_H
#define MATERN_COV_H

#include <Eigen/Dense>
#include <vector>

// Function to compute the Matern covariance matrix
Eigen::MatrixXd matern_cov_spat(double sigmasq, double range_, double v, const Eigen::MatrixXd& input_coords);

// Negative log-likelihood function
double neg_log_likelihood_nugget(const Eigen::MatrixXd& input_coords, const Eigen::VectorXd& y, double params);

#endif  // MATERN_COV_H#pragma once
