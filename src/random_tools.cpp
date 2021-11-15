#include "random_tools.h"

namespace bgv
{
std::default_random_engine RandomGenerator::generator = std::default_random_engine();

Eigen::MatrixXi RandomGenerator::generateUniformIntegerMatrix(int n_rows, int n_cols, int min, int max)
{
    std::uniform_int_distribution<int> distribution(min, max);
    
    Eigen::MatrixXi m = Eigen::MatrixXi::NullaryExpr(n_rows, n_cols, [&](){return distribution(generator);});
    return m;
}


Eigen::MatrixXi RandomGenerator::generateErrorMatrix(int n_rows, int n_cols, int min, int max)
{
    //should be gaussian (or rather binomial ???) here, will use uniform distribution for the time being
    std::uniform_int_distribution<int> distribution(std::max(-16, min), std::min(16, max));
    
    Eigen::MatrixXi m = Eigen::MatrixXi::NullaryExpr(n_rows, n_cols, [&](){return distribution(generator);});
    return m;
}

}