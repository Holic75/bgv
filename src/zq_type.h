#pragma once
#include <iostream>
#include <assert.h>

namespace bgv
{
//assume that q is even
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


    static int restrict(int z, int modulus)
    {
        //restrict value to [-modulus/2;  modulus/2]
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
};
}

