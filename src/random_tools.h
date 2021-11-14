#pragma once
#include <random>
#include <Eigen/Dense>

namespace bgv
{

class RandomGenerator
{
    static std::default_random_engine generator;
public:
    static int generateUniformInteger(int min, int max) 
    {
        std::uniform_int_distribution<int> distribution(min, max);
        return distribution(generator);
    }

    static Eigen::MatrixXi generateUniformIntegerMatrix(int n_rows, int n_cols, int min, int max)
    {
        std::uniform_int_distribution<int> distribution(min, max);
        
        Eigen::MatrixXi m = Eigen::MatrixXi::NullaryExpr(n_rows, n_cols, [&](){return distribution(generator);});
        return m;
    }


    static Eigen::MatrixXi generateErrorMatrix(int n_rows, int n_cols, int min, int max)
    {
        //should be gaussian here, will use uniform distribution for the time being
        std::uniform_int_distribution<int> distribution(-16, 16);
        
        Eigen::MatrixXi m = Eigen::MatrixXi::NullaryExpr(n_rows, n_cols, [&](){return distribution(generator);});
        return m;
    }
};

}




