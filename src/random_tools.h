#pragma once
#include <random>
#include <Eigen/Dense>

namespace bgv
{

/// \brief Some utility functions for random number generation
/// 
/// Always uses the same default seed to ensure repeatability of tests
/// TODO: add option to provide a specific seed
class RandomGenerator
{
    static std::default_random_engine generator;
public:

    /// \brief Generate random integer from uniform distribution on [min, max]
	///
	/// @param min Left bound.
	/// @param max Right Bound.
	/// @return Random integer from uniform distribution on [min, max].
    static int generateUniformInteger(int min, int max) 
    {
        std::uniform_int_distribution<int> distribution(min, max);
        return distribution(generator);
    }


    /// \brief Generate random matrix of integers from uniform distribution on [min, max]
	///
    /// @param n_rows Number of rows in matrix.
	/// @param n_cols Number of columns in matrix.
	/// @param min Left bound.
	/// @param max Right Bound.
	/// @return Random matrix from uniform distribution on [min, max].
    static Eigen::MatrixXi generateUniformIntegerMatrix(int n_rows, int n_cols, int min, int max);

    /// \brief Generate random erro matrix
	///
    /// @param n_rows Number of rows in matrix.
	/// @param n_cols Number of columns in matrix.
	/// @param min Left bound on error value.
	/// @param max Right Bound on error value.
	/// @return Random matrix from uniform distribution on [min, max].
    /// TODO: replace uniform distribution with proper one
    static Eigen::MatrixXi generateErrorMatrix(int n_rows, int n_cols, int min, int max);
};

}




