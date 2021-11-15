#include "encryption.h"

namespace bgv
{

Eigen::MatrixXi EncryptionEngine::encrypt(const Eigen::VectorXi& s, int modulus, const Eigen::VectorXi& m)
{
    Eigen::MatrixXi a(m.size(), s.size() + 1);
    a.block(0, 1, m.size(), s.size()) = RandomGenerator::generateUniformIntegerMatrix(m.size(), s.size(), -modulus/2, modulus/2);

    Eigen::MatrixXi e = RandomGenerator::generateErrorMatrix(m.size(), 1, -modulus/4, modulus/4);
    Eigen::VectorXi b = a.block(0, 1, m.size(), s.size()) * s + 2 * e + m;
    b = ZqElement::restrict(b, modulus);
    a.col(0) = b;
    return a;
}


Eigen::MatrixXi EncryptionEngine::encrypt(const Eigen::VectorXi& s, int modulus, int value)
{
    Eigen::VectorXi m(1);
    m(0) = value;
    return encrypt(s, modulus, m);
}


Eigen::MatrixXi EncryptionEngine::encryptKeyPowersOf2(const Eigen::VectorXi& new_s, int modulus, const Eigen::VectorXi& old_s) 
{
    //here we want to encrypt 2^tau * (1, s[i]) x (1, s[j])
    //we will need a total of log(q) * (1 + s.size()) * (2 + s.size()) / 2 encryptions
    //so random matrix will need to have size: {log(q) * (1 + old_s.size())} * new_s.size() 
    int log2q = log2(modulus) - 1;
    Eigen::VectorXi s1(old_s.size() + 1);
    s1(0) = 1;
    s1.tail(old_s.size()) = old_s;
    
    Eigen::MatrixXi ss_mat = s1 * s1.transpose();
    Eigen::MatrixXi ss_vec = storeSymmetricMatrixAsVector(ss_mat, false);

    auto ss_powers_of2_vec = createPowersOf2Vector(ss_vec, log2q);

    return encrypt(new_s, modulus, ss_powers_of2_vec);
}


Eigen::MatrixXi EncryptionEngine::encryptKeyBitsPowersOf2(const Eigen::VectorXi& new_s, int new_modulus, const Eigen::VectorXi& old_s, int old_modulus) 
{
    assert(new_modulus < old_modulus);
    assert(new_modulus > 0);
    //here we want to encrypt 2^k * {(1, s[i]) x (1, s[j])}_t
    //where (1, s[i])(1, s[j]) = sum_t 2^t (1, s[i]) x (1, s[j])}_t
    //we will need a total of log(q_new) * log(q_old) * (1 + s.size()) * (2 + s.size()) / 2 encryptions

    int log2q_old = log2(old_modulus) - 1;
    int log2q_new = log2(new_modulus) - 1;

    Eigen::VectorXi s1(old_s.size() + 1);
    s1(0) = 1;
    s1.tail(old_s.size()) = old_s;
    Eigen::MatrixXi ss_mat = s1 * s1.transpose();
    Eigen::VectorXi ss_vec = storeSymmetricMatrixAsVector(ss_mat, false);
    ss_vec = createVectorBitDecomposition(ss_vec, log2q_old);
    auto ss_powers_of2_vec = createPowersOf2Vector(ss_vec, log2q_new);

    return encrypt(new_s, new_modulus, ss_powers_of2_vec);
}


Eigen::MatrixXi EncryptionEngine::refresh(const Eigen::VectorXi& h, const Eigen::MatrixXi& info, int old_modulus, int new_modulus)
{
    assert(new_modulus < old_modulus);
    assert(new_modulus > 0);
    //take c = sum_ij h_ij * s[i]s[j] 
    //= sum_ijt h_ij 2^t {s[i]s[j]}_t   (i.e. decomposition with respect to binary representation of the key)
    //apply modulus switching
    //u_ijt = {h_ij 2^t} -> v_ijt = (p/q) * u_ijt s.t.  u_ijt mod 2 = v_ijt mod 2
    //c = sum_ijt v_ijt {s[i]s[j]}_t = sum_ijtk v_ijtk * 2^k {s[i]s[j]}_t
    int log2_q_old = log2(old_modulus) - 1;
    int half_new_modulus  = new_modulus / 2;
    Eigen::VectorXi h_powers_of2 = ZqElement::restrict(createPowersOf2Vector(h, log2_q_old), old_modulus);

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

    Eigen::MatrixXi c = relinearize(h_powers_of2, info);
    return c;
}


Eigen::MatrixXi EncryptionEngine::multiplyCyphertextsLong(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2, int modulus)
{
    //c1 and c2 are row vectors
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

    h = ZqElement::restrict(h, modulus);
    Eigen::VectorXi h_vec = ZqElement::restrict(storeSymmetricMatrixAsVector(h, true), modulus);

    return h_vec;
}


Eigen::MatrixXi EncryptionEngine::multiplyCyphertexts(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2, Eigen::MatrixXi info, int modulus)
{
    auto h_vec =  multiplyCyphertextsLong(c1, c2, modulus);
    Eigen::MatrixXi c = relinearize(h_vec, info);

    return c;
}


Eigen::MatrixXi EncryptionEngine::relinearize(const Eigen::VectorXi& h_vec,  Eigen::MatrixXi info)
{
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


Eigen::VectorXi EncryptionEngine::storeSymmetricMatrixAsVector(const Eigen::MatrixXi& m, bool sum)
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

    return vec;
}


Eigen::MatrixXi EncryptionEngine::addCyphertexts(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2, int modulus)
{
    //c1 and c2 are row vectors
    assert(c1.size() == c2.size());
    assert(c1.rows() == 1);

    Eigen::MatrixXi c_sum = c1 + c2;
    return ZqElement::restrict(c_sum, modulus);
}


Eigen::VectorXi EncryptionEngine::addCyphertextsLong(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2, int modulus)
{
    //c1 and c2 are row vectors
    //c_i = (<a_i, s>, a_i), so size(c_i) = size(s) + 1
    //m_1 + m_2 = (b_1 - <a_1, s>) + (b_2 - <a_2,s>) = b1 + b2 - <a_1 + a_2,s>
    //= sum_i h_i * s_i

    assert(c1.size() == c2.size());
    assert(c1.rows() == 1);
    Eigen::VectorXi c_sum = addCyphertexts(c1, c2, modulus).transpose();
    c_sum.tail(c1.cols() - 1) = -c_sum.tail(c1.cols() - 1);

    Eigen::VectorXi h_vec = Eigen::VectorXi::Zero((c1.cols() + 1) * c1.cols() / 2);
    h_vec.head(c1.cols()) = c_sum;
    return h_vec;
}


Eigen::VectorXi EncryptionEngine::decrypt(const Eigen::VectorXi& s, int modulus, const Eigen::MatrixXi& cyphertext)
{
    Eigen::VectorXi result = (cyphertext.block(0, 0, cyphertext.rows(), 1) 
                        - cyphertext.block(0, 1, cyphertext.rows(), s.size()) * s);
    result = ZqElement::restrict(result , modulus);
    result = result.unaryExpr([](int elem){return abs(elem % 2);});
    return result;
}


Eigen::VectorXi EncryptionEngine::createVectorBitDecomposition(const Eigen::VectorXi& v, int n_bits)
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


Eigen::VectorXi EncryptionEngine::createPowersOf2Vector(const Eigen::MatrixXi& vec, int n_bits)
{
    auto powers_of2 = createPowersOf2(n_bits);
    
    Eigen::MatrixXi powers_of2_mat = vec * powers_of2.transpose();
    Eigen::VectorXi powers_of2_vec 
            = Eigen::Map<Eigen::VectorXi>(powers_of2_mat.data(), powers_of2_mat.cols() * powers_of2_mat.rows());
    
    return powers_of2_vec;
}


Eigen::VectorXi EncryptionEngine::createPowersOf2(int size)
{
    Eigen::VectorXi v(size);
    v(0) = 1;
    for (int i = 1; i < size; i++)
    {
        v(i) = (v(i - 1) << 1);
    }

    return  v;
}

}