#pragma once
#include "random_tools.h"
#include "zq_type.h"
#include <Eigen/Dense>

namespace bgv
{

class EncryptionEngine
{

public:
    static Eigen::MatrixXi encrypt(const Eigen::VectorXi& s, int modulus, bool value)
    {
        Eigen::MatrixXi a(1, s.size() + 1);
        a.block(0, 1, 1, s.size()) = RandomGenerator::generateUniformIntegerMatrix(1, s.size(), -modulus/2, modulus/2);
        
        Eigen::MatrixXi e = RandomGenerator::generateErrorMatrix(1, 1, -modulus/4, modulus/4);

        int b = (a.block(0, 1, 1, s.size()) * s + 2 * e)(0, 0);
        b = ZqElement::restrict(b + (value ? 1 : 0), modulus);
        a(0, 0) = b;
        return a; //will return a row vector
    }


    static Eigen::MatrixXi encryptKeyPowersOf2(const Eigen::VectorXi& new_s, int modulus, const Eigen::VectorXi& old_s) 
    {
        //here we want to encrypt 2^tau * (1, s[i]) x (1, s[j])
        //we will need a total of log(q) * (1 + s.size()) * (2 + s.size()) / 2 encryptions
        //so random matrix will need to have size: {log(q) * (1 + old_s.size())} * new_s.size() 
        int log2q = log2(modulus) - 1;
        //int encryption_size = log2q * (1 + old_s.size()) * (1 + old_s.size()); //* (2 + old_s.size()) / 2;
        int encryption_size = log2q * (1 + old_s.size()) * (2 + old_s.size()) / 2;

        Eigen::MatrixXi a(encryption_size, new_s.size() + 1);
        a.block(0, 1, encryption_size, new_s.size()) = RandomGenerator::generateUniformIntegerMatrix(encryption_size, new_s.size(), -modulus/2, modulus/2);

        Eigen::MatrixXi e = RandomGenerator::generateErrorMatrix(encryption_size, 1, -modulus/4, modulus/4);

        Eigen::VectorXi s1(old_s.size() + 1);
        s1(0) = 1;
        s1.tail(old_s.size()) = old_s;
        
        Eigen::MatrixXi ss_mat = s1 * s1.transpose();
        Eigen::MatrixXi ss_vec = storeSymmetricMatrixAsVector(ss_mat, false);

        auto ss_powers_of2_vec = createPowersOf2Vector(ss_vec, log2q);

        Eigen::VectorXi b = (a.block(0, 1, encryption_size, new_s.size()) * new_s + 2 * e);  

        b = (b + (ss_powers_of2_vec)).unaryExpr([modulus](int elem){return ZqElement::restrict(elem, modulus);});
        a.block(0, 0, encryption_size, 1) = b;
        return a;
    }


    static Eigen::MatrixXi encryptKeyBitsPowersOf2(const Eigen::VectorXi& new_s, int modulus, const Eigen::VectorXi& old_s) 
    {
        //here we want to encrypt 2^k * {(1, s[i]) x (1, s[j])}_t
        //where (1, s[i])(1, s[j]) = sum_t 2^t (1, s[i]) x (1, s[j])}_t
        //we will need a total of log(q) * log(q) * (1 + s.size()) * (2 + s.size()) / 2 encryptions

        int log2q = log2(modulus) - 1;
        int encryption_size = log2q * log2q * (1 + old_s.size()) * (2 + old_s.size()) / 2;

        Eigen::MatrixXi a(encryption_size, new_s.size() + 1);
        a.block(0, 1, encryption_size, new_s.size()) = RandomGenerator::generateUniformIntegerMatrix(encryption_size, new_s.size(), -modulus/2, modulus/2);

        Eigen::MatrixXi e = RandomGenerator::generateErrorMatrix(encryption_size, 1, -modulus/4, modulus/4);

        Eigen::VectorXi s1(old_s.size() + 1);
        s1(0) = 1;
        s1.tail(old_s.size()) = old_s;
        
        Eigen::MatrixXi ss_mat = s1 * s1.transpose();
        Eigen::VectorXi ss_vec = storeSymmetricMatrixAsVector(ss_mat, false);
        ss_vec = createVectorBitDecomposition(ss_vec, log2q);
        auto ss_powers_of2_vec = createPowersOf2Vector(ss_vec, log2q);

        Eigen::VectorXi b = (a.block(0, 1, encryption_size, new_s.size()) * new_s + 2 * e);  

        b = (b + (ss_powers_of2_vec)).unaryExpr([modulus](int elem){return ZqElement::restrict(elem, modulus);});
        a.block(0, 0, encryption_size, 1) = b;
        return a;
    }

    static Eigen::VectorXi refresh(const Eigen::VectorXi& h, const Eigen::MatrixXi& info, int old_modulus, int new_modulus)
    {
        assert(new_modulus < old_modulus);
        assert(new_modulus > 0);
        //take c = sum_ij h_ij * s[i]s[j] 
        //= sum_ijt h_ij 2^t {s[i]s[j]}_t   (i.e. decomposition with respect to binary representation of the key)
        //apply modulus switching
        //u_ijt = {h_ij 2^t} -> v_ijt = (p/q) * u_ijt s.t.  u_ijt mod 2 = v_ijt mod 2
        //c = sum_ijt v_ijt {s[i]s[j]}_t = sum_ijtk v_ijtk * 2^k {s[i]s[j]}_t
        int log2_q_old = log2(old_modulus);
        int half_new_modulus  = new_modulus / 2;
        Eigen::VectorXi h_powers_of2 = createPowersOf2Vector(h, log2_q_old).unaryExpr([old_modulus](int elem){return ZqElement::restrict(elem, old_modulus);});

        for (int i = 0; i < h_powers_of2.size(); i++)
        {
            int old_mod = h_powers_of2(i) % 2;
            int new_value = (h_powers_of2(i) * new_modulus / old_modulus);
            int new_mod = new_value % 2;
            // new_mod and old_mod will have same signs
            if (old_mod > new_mod)
            {
                new_value = (new_mod == half_new_modulus) ? new_value - 1 : new_value + 1;
            }
            else if (old_mod < new_mod)
            {
                new_value = (new_mod == - half_new_modulus) ? new_value + 1 : new_value - 1;
            }
            h_powers_of2(i) = new_value;
        }

        int log2q_new = log2(new_modulus);
        auto h_vec_binary = createVectorBitDecomposition(h_powers_of2, log2q_new);
        Eigen::MatrixXi c = h_vec_binary.transpose() * info;
        return c;
    }



    static Eigen::MatrixXi multiplyCyphertexts(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2, Eigen::MatrixXi info, int modulus)
    {
        //c1 and c2 are row vectoes
        assert(c1.size() == c2.size());
        assert(c1.rows() == 1);
        //c1 = (b,  a ) = (<a ,s> + 2e +  m,  a )
        //c2 = (b', a') = (<a',s> + 2e' + m', a')
        //m*m'  ~ (b - <a,s>) * (b' - <a',s>) =  sum_ij h_ij (1, s)[i]*(1, s)[j] = sum_ijt h_ijt * {2^t (1, s)[i] * (1, s)[j]}
        //info = {2^t (1, s)[i] * (1, s)[j]}
        //h_i,j = b * b', if i = j = 0
        //      = -b * a'[j - 1],  if i = 0, j > 0
        //      = a[i - 1] * b', if i > 0, j = 0
        //      = a[i - 1] * a'[j - 1], if i > 0, j > 0
        
        Eigen::MatrixXi h(c1.cols(), c1.cols());
        h(0, 0) = c1(0, 0) * c2(0, 0);
        h.block(1, 0, c1.cols() - 1, 1) = -c2(0, 0) * c1.block(0, 1, 1, c2.cols() - 1).transpose();
        h.block(0, 1, 1, c1.cols() - 1) = -c1(0, 0) * c2.block(0, 1, 1, c1.cols() - 1);
        h.block(1, 1, c1.cols() - 1, c1.cols() - 1) = c1.block(0, 1, 1, c2.cols() - 1).transpose() * c2.block(0, 1, 1, c2.cols() - 1);

        h = h.unaryExpr([modulus](int elem){return ZqElement::restrict(elem, modulus);});
        Eigen::MatrixXi h_vec = storeSymmetricMatrixAsVector(h, true).unaryExpr([modulus](int elem){return ZqElement::restrict(elem, modulus);});

        //auto h_vec =  Eigen::Map<Eigen::VectorXi>(h.data(), h.cols() * h.rows());
        int n = h_vec.size();
        int log2q = info.rows() / n;
        auto h_vec_binary = createVectorBitDecomposition(h_vec, log2q);

        //sum_ij h_ij (1, s)[i]*(1, s)[j] = sum_ijt hb_ijt * {2^t (1, s)[i] * (1, s)[j]} = 
        //sum_ijt hb_ijt * {b_ijt - <a_ijt, new_s>} = sum_ijt hb_ijt * b_ijt-  <sum_ijt hb_ijt * a_ijt, new_s>
        //now we can aquire cyphertext encryption in new key as (hb_ijt * b_ijt, sum_ijt hb_ijt * a_ijt)
        //it seems that in this case the error should be bounded by (old_s.size() + 1)^2 * (log_2(q) - 1)
        //if we do not use binary decomposition then it is not bounded at all (since h_ij can normally take any value in Z_q)
        Eigen::MatrixXi c = h_vec_binary.transpose() * info;
        return c;
    }


    static Eigen::VectorXi storeSymmetricMatrixAsVector(const Eigen::MatrixXi& m, bool sum)
    {
        assert(m.rows() == m.cols());
        Eigen::VectorXi vec(m.cols() * (m.cols() + 1) / 2);
        int idx_start = 0;
        for (int i = 0; i < m.cols(); i++)
        {
            int elts_in_col = m.rows() - i;
            vec.segment(idx_start + 1, elts_in_col - 1) = m.col(i).segment(i + 1, elts_in_col - 1);
            if (sum)
            {
                vec.segment(idx_start + 1, elts_in_col - 1) += m.row(i).segment(i + 1, elts_in_col - 1).transpose();
            }
            vec(idx_start) = m(i, i);
            idx_start += elts_in_col;
        }
        assert(idx_start == vec.size());
        //auto mm = m;
        //vec = Eigen::Map<Eigen::VectorXi>(mm.data(), mm.cols() * mm.rows());

        return vec;
    }


    static Eigen::MatrixXi addCyphertexts(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2)
    {
        //c1 and c2 are row vectors
        assert(c1.size() == c2.size());
        assert(c1.rows() == 1);

        return c1 + c2;
    }


    static Eigen::VectorXi createVectorBitDecomposition(const Eigen::VectorXi& v, int n_bits)
    {
        Eigen::VectorXi vb(n_bits * v.size());
        Eigen::VectorXi vv = v;

        for (int i = 0; i < n_bits; i++)
        {
            vb.segment(i* v.size(), v.size()) = vv.unaryExpr([](int elem)
                                                             {
                                                                 if (elem > 0)
                                                                 {
                                                                     return elem & 1;
                                                                 }
                                                                 else
                                                                 {
                                                                     return -(abs(elem) & 1);
                                                                 }
                                                             });
            vv = vv.unaryExpr([](int elem)
                            {
                                if (elem > 0)
                                {
                                    return elem >> 1;
                                }
                                else
                                {
                                    return -(abs(elem) >> 1);
                                }
                            });
        }
        
        return vb;
    }


    static bool decrypt(const Eigen::VectorXi& s, int modulus, const Eigen::MatrixXi& cyphertext)
    {
        int result = (cyphertext.block(0, 0, cyphertext.rows(), 1) 
                            - cyphertext.block(0, 1, cyphertext.rows(), s.size()) * s)(0, 0);
        result = ZqElement::restrict(result , modulus);
        result = result % 2;
        return result;
    }


    static Eigen::VectorXi createPowersOf2Vector(const Eigen::MatrixXi& vec, int n_bits)
    {
        auto powers_of2 = createPowersOf2(n_bits);
        
        Eigen::MatrixXi powers_of2_mat = vec * powers_of2.transpose();
        Eigen::VectorXi powers_of2_vec 
                = Eigen::Map<Eigen::VectorXi>(powers_of2_mat.data(), powers_of2_mat.cols() * powers_of2_mat.rows());
        
        return powers_of2_vec;
    }



    static  Eigen::VectorXi createPowersOf2(int size)
    {
        Eigen::VectorXi v(size);
        v(0) = 1;
        for (int i = 1; i < size; i++)
        {
            v(i) = (v(i - 1) << 1);
        }

        return  v;
    }


    static int log2(unsigned int q)
    {
        int count = 0;
        while (q > 0)
        {
            count ++;
            q = q >> 1;
        }
        return count;
    }
    //addition of cypertexts is trivial
    //multiplication will require extra information about s x s decomposition
};

}