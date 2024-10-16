#include "matern_cov.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions> // For advanced matrix operations
#include <gsl/gsl_sf_bessel.h>               // GSL library for special functions (e.g., Bessel)
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846  // Pi constant to a high precision
#endif

// Function to compute the Matern covariance matrix
Eigen::MatrixXd matern_cov_spat(double sigmasq, double range_, double v, const Eigen::MatrixXd& input_coords) {
    int n = input_coords.rows();
    Eigen::MatrixXd cov_matrix(n, n);

	#pragma omp parallel for
    // Calculate pairwise distances
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double dist = (input_coords.row(i) - input_coords.row(j)).norm();
            if (dist == 0) {
                cov_matrix(i, j) = sigmasq;
            }
            else {
                // Mat?rn covariance formula using Bessel function
                double factor = sigmasq * std::pow(2, 1 - v) / tgamma(v);
                double cov = factor * std::pow(dist / range_, v) * gsl_sf_bessel_Knu(v, dist / range_);
                cov_matrix(i, j) = cov;
                cov_matrix(j, i) = cov;  // Symmetric matrix
            }
        }
    }

    return cov_matrix;
}

// Function to calculate the negative log-likelihood
double neg_log_likelihood_nugget(const Eigen::MatrixXd& input_coords, const Eigen::VectorXd& y, double params) {
    // Compute covariance matrix
    Eigen::MatrixXd cov_matrix = matern_cov_spat(1.0, 2.0, 0.5, input_coords);

    // Add a jitter term (nugget) for numerical stability
    cov_matrix += Eigen::MatrixXd::Identity(cov_matrix.rows(), cov_matrix.cols()) * params;

    // Cholesky decomposition
    Eigen::LLT<Eigen::MatrixXd> llt(cov_matrix);
    Eigen::MatrixXd L = llt.matrixL();

    // Log determinant
    double log_det = 2 * L.diagonal().array().log().sum();

    // Calculate beta
    Eigen::VectorXd beta = (input_coords.transpose() * llt.solve(input_coords)).inverse()
        * input_coords.transpose() * llt.solve(y);
    // Mean
    Eigen::VectorXd mu = input_coords * beta;
    Eigen::VectorXd y_mu = y - mu;

    Eigen::VectorXd alpha = L.triangularView<Eigen::Lower>().solve(y_mu);
    double quad_form = alpha.transpose() * alpha;

    // Negative log-likelihood
    // int n = y.size();
    int n = static_cast<int>(input_coords.rows());   //  Added explicit type conversions using static_cast<int>

    double neg_log_lik = 0.5 * (n * std::log(2 * M_PI) + log_det + quad_form);

    return neg_log_lik;
}