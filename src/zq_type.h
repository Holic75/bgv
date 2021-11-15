#pragma once
#include <iostream>
#include <assert.h>
#include <Eigen/Dense>

namespace bgv
{
/// \brief Utility class implementing Zq elements restricted to [-q/2, q/2]
/// 
/// q is assumed to be odd
class ZqElement
{
private:
    int value_;
    const int q_;
    void restrict_()
    {
        value_ = restrict(value_, q_);
    }

public:
    ZqElement(int q, int value = 0)
        :value_(value), q_(q)
    {
        assert(q % 2 == 1 && q > 0);
        restrict_();
    };
    operator int() const { return value_; }
    ZqElement operator+ (ZqElement other) {return ZqElement(this->q_, this->value_ + other.value_);}
    ZqElement operator- (ZqElement other) {return ZqElement(this->q_, this->value_ - other.value_);}
    ZqElement operator* (ZqElement other) {return ZqElement(this->q_, this->value_ * other.value_);}
    ZqElement operator- () {return ZqElement(this->q_, -this->value_);};
    friend std::ostream & operator << (std::ostream &out, const ZqElement z)
    {
        return out << z.value_;
    }

    /// \brief Restrict number into [-modulus/2, modulus/2] interval.
	///
    /// @param z Number.
	/// @param modulus Modulus.
	/// @return z restricted to [-modulus/2, modulus/2] interval
    static int restrict(int z, int modulus)
    {
        z = z % modulus;
        if (z > modulus/2)
        {
            z -= modulus;
        }
        else if (z < - modulus/2)
        {
            z += modulus;
        }
        return z;
    }

    /// \brief Restrict matrix elements into [-modulus/2, modulus/2] interval.
	///
    /// @param z Matrix.
	/// @param modulus Modulus.
	/// @return Matrix with elements restricted to [-modulus/2, modulus/2] interval
    template<int R, int C>
    static Eigen::Matrix<int, R, C> restrict(const Eigen::Matrix<int, R, C>& z, int modulus)
    {
        return z.unaryExpr([modulus](int elem){return restrict(elem, modulus);});
    }
};
}

