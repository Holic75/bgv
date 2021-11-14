#pragma once
#include "random_tools.h"
#include "zq_type.h"
#include <Eigen/Dense>

namespace bgv
{

class EncryptionEngine
{

public:
    static Eigen::MatrixXi encrypt(const Eigen::VectorXi& s, int modulo, bool value)
    {
        Eigen::MatrixXi a(1, s.size() + 1);
        a.block(0, 1, 1, s.size()) = RandomGenerator::generateUniformIntegerMatrix(1, s.size(), -modulo/2, modulo/2);
        
        Eigen::MatrixXi e = RandomGenerator::generateErrorMatrix(1, 1, -modulo/4, modulo/4);

        int b = (a.block(0, 1, 1, s.size()) * s + 2 * e)(0, 0);
        b = ZqElement::restrict(b + (value ? 1 : 0), modulo);
        a(0, 0) = b;
        return a; //will return a row vector
    }


    static Eigen::MatrixXi encryptKeyPowersOf2(const Eigen::VectorXi& new_s, int modulo, const Eigen::VectorXi& old_s) 
    {
        //here we want to encrypt 2^tau * (1, s[i]) x (1, s[j])
        //we will need a total of log(q) * (1 + s.size()) * (2 + s.size()) / 2 encryptions
        //so random matrix will need to have size: {log(q) * (1 + old_s.size())} * new_s.size() 
        int log2q = log2(modulo) - 1;
        //int encryption_size = log2q * (1 + old_s.size()) * (1 + old_s.size()); //* (2 + old_s.size()) / 2;
        int encryption_size = log2q * (1 + old_s.size()) * (2 + old_s.size()) / 2;

        Eigen::MatrixXi a(encryption_size, new_s.size() + 1);
        a.block(0, 1, encryption_size, new_s.size()) = RandomGenerator::generateUniformIntegerMatrix(encryption_size, new_s.size(), -modulo/2, modulo/2);

        Eigen::MatrixXi e = RandomGenerator::generateErrorMatrix(encryption_size, 1, -modulo/4, modulo/4);

        Eigen::VectorXi s1(old_s.size() + 1);
        s1(0) = 1;
        s1.tail(old_s.size()) = old_s;
        
        Eigen::MatrixXi ss_mat = s1 * s1.transpose();
        Eigen::MatrixXi ss_vec = storeSymmetricMatrixAsVector(ss_mat, false);

        auto powers_of2_vec = createPowersOf2Vector(log2q);
        
        Eigen::MatrixXi ss_powers_of2_mat = ss_vec* powers_of2_vec.transpose();
        Eigen::VectorXi ss_powers_of2_vec 
                = Eigen::Map<Eigen::VectorXi>(ss_powers_of2_mat.data(), ss_powers_of2_mat.cols() * ss_powers_of2_mat.rows());

        Eigen::VectorXi b = (a.block(0, 1, encryption_size, new_s.size()) * new_s + 2 * e);  

        b = (b + (ss_powers_of2_vec)).unaryExpr([modulo](int elem){return ZqElement::restrict(elem, modulo);});
        a.block(0, 0, encryption_size, 1) = b;
        return a;
    }



    static Eigen::MatrixXi multiplyCyphertexts(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2, Eigen::MatrixXi info, int modulo)
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

        h = h.unaryExpr([modulo](int elem){return ZqElement::restrict(elem, modulo);});
        Eigen::MatrixXi h_vec = storeSymmetricMatrixAsVector(h, true).unaryExpr([modulo](int elem){return ZqElement::restrict(elem, modulo);});

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


    static bool decrypt(const Eigen::VectorXi& s, int modulo, const Eigen::MatrixXi& cyphertext)
    {
        int result = (cyphertext.block(0, 0, cyphertext.rows(), 1) 
                            - cyphertext.block(0, 1, cyphertext.rows(), s.size()) * s)(0, 0);
        result = ZqElement::restrict(result , modulo);
        result = result % 2;
        return result;
    }


    static  Eigen::VectorXi createPowersOf2Vector(int size)
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