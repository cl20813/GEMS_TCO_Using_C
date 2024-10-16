#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <unordered_map>
#include <iomanip> // For std::get_time

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <unsupported/Eigen/MatrixFunctions> // For advanced matrix operations
#include "matern_cov.h"
#include <omp.h>


// Define a custom structure
struct DataFrame {
    Eigen::MatrixXd data;
    std::vector<std::string> orbits;
};

// Function to read the CSV data and return the custom DataFrame
DataFrame read_csv(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;

    // Get the number of rows and columns
    int rows = 0;
    int cols = 0;
    std::getline(file, line);  // Read the header
    std::stringstream header_stream(line);
    std::string cell;
    while (std::getline(header_stream, cell, ',')) {
        cols++;
    }

    while (std::getline(file, line)) {
        rows++;
    }

    // Reset file stream
    file.clear();
    file.seekg(0, std::ios::beg);
    std::getline(file, line); // Skip header

    // Initialize DataFrame
    DataFrame df;
    df.data.resize(rows, cols - 1);  // Exclude the 'Orbit' column initially

    int row = 0;
    while (std::getline(file, line)) {
        std::stringstream line_stream(line);
        std::vector<std::string> row_cells;
        while (std::getline(line_stream, cell, ',')) {
            row_cells.push_back(cell);
        }

        std::string orbit;
        for (int col = 0; col < cols; ++col) {
            if (col == 3) {  // Assuming 'Orbit' is the fourth column
                // Parse the date-time value using std::get_time
                std::istringstream ss(row_cells[col]);
                std::tm tm = {};
                ss >> std::get_time(&tm, "%Y-%m-%d  %H:%M");
                if (ss.fail()) {
                    std::cerr << "Failed to parse date-time: " << row_cells[col] << std::endl;
                }
                else {
                    std::ostringstream oss;
                    oss << std::put_time(&tm, "%Y-%m-%d  %H:%M");
                    orbit = oss.str();
                }
            }
            else {
                df.data(row, col < 3 ? col : col - 1) = std::stod(row_cells[col]);
            }
        }
        df.orbits.push_back(orbit);
        row++;
    }
    return df;
}

// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

};

struct my_functor : Functor<double>
{
    const Eigen::MatrixXd& input_coords;  // Coordinates (latitude and longitude) so nx2 matrix
    const Eigen::VectorXd& y;             // Response variable (ColumnAmountO3)  nx1 vector

    // Initialize an object 
    my_functor(const Eigen::MatrixXd& input_coords_, const Eigen::VectorXd& y_) // constructor of my_functor class with two arguments
        : Functor<double>(1, 1), input_coords(input_coords_), y(y_) {}  //   the functor has 1 input and 1 output, and then two references are initialized


    // Operator to compute the residuals (the negative log-likelihood)
    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
	{
    fvec(0) = neg_log_likelihood_nugget(input_coords, y, x(0));  // Nugget parameter
    return 0;
	}

};


int main() {
    std::string file_path = "C:\\Users\\joonw\\TCO\\data_engineering\\sparse_data_24_7_300_N3035_E100110.csv";

    // Read the data into a custom DataFrame
    DataFrame df = read_csv(file_path);

    // Ensure the data is read correctly
    // std::cout << "Rows read: " << df.data.rows() << std::endl;
    // std::cout << "Columns read: " << df.data.cols() << std::endl;
    // std::cout << "Orbit values read: " << df.orbits.size() << std::endl;

    // Extract unique keys for 'Orbit'
    std::vector<std::string> unique_orbits = df.orbits;           // without this I am doing in-place
    std::sort(unique_orbits.begin(), unique_orbits.end());
    unique_orbits.erase(std::unique(unique_orbits.begin(), unique_orbits.end()), unique_orbits.end());

    std::vector<int> results(unique_orbits.size());

    // Map to organize data by orbit
    std::unordered_map<std::string, std::vector<Eigen::VectorXd>> orbit_map;


    for (int i = 0; i < df.data.rows(); ++i) {
        orbit_map[df.orbits[i]].push_back(df.data.row(i));
    }

    // Perform your optimization for each orbit

	#pragma omp parallel for
    for (const std::string& orbit_key : unique_orbits) { // Using unique orbits here
        const std::vector<Eigen::VectorXd>& cur_data_vec = orbit_map[orbit_key];

        // Convert cur_data_vec to Eigen::MatrixXd
        Eigen::MatrixXd cur_data(cur_data_vec.size(), df.data.cols());
        for (size_t j = 0; j < cur_data_vec.size(); ++j) {
            cur_data.row(j) = cur_data_vec[j];
        }

        // Separate inputs and outputs
        Eigen::MatrixXd input_coords = cur_data.leftCols(2); // Assuming the first two columns are coordinates
        Eigen::VectorXd y = cur_data.col(2); // Assuming the third column is the target variable

        // std::cout << "input_coords: " << input_coords << std::endl;
        // std::cout << "y: " << y << std::endl;

        // Initial guess for the parameter
        Eigen::VectorXd params(1);  // initializes an Eigen vector with size 1
        params << 0.1; // Initial value for the parameter

        // Optimization using Levenberg-Marquardt algorithm

        my_functor functor(input_coords, y);

        // Numerical differentiation for approximating the Jacobian
        Eigen::NumericalDiff<my_functor> numDiff(functor);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<my_functor>, double> lm(numDiff);
  
        lm.parameters.maxfev = 200;  // Reduced from 2000 for efficiency
        lm.parameters.xtol = 1.0e-5;  // Looser tolerance for faster convergence 1.0e-7;  // or 1.0e-10

        std::cout << lm.parameters.maxfev << std::endl;

        // Perform the optimization
        int ret = lm.minimize(params);

        std::cout << "number of iterations: "<< lm.iter << std::endl;  
        std::cout << "return status: " << ret << std::endl;  // 1: Successful convergence. 0 : Hit the maximum number of iterations without converging.- 1 : Encountered an error or failed to converge.

        std::cout << "Optimized parameter for Orbit " << orbit_key << "is " << params(0) << std::endl;

        // std::cout << "Shape of cur_data for Orbit " << orbit_key << ": " << cur_data.rows() << " rows" << std::endl;

    }

    return 0;
}


