#include "zq_type.h"
#include "random_tools.h"
#include "encryption.h"
#include <Eigen/Dense>

int main (int argc, char *argv[]) 
{
    /*bgv::ZqElement z1(11, 5), z2(11, 3);

    Eigen::MatrixXi m(2,2);
    m(0, 0) =  z2 - z1;
    m(1, 0) = z1 + z2;
    m(0, 1) = z1 * z2;
    m(1, 1) = -z1;*/


    /*Eigen::VectorXi s1_test(2);
    s1_test << -1, 2;

    Eigen::VectorXi s2_test(2);
    s2_test << 0, -1;

    int q_test = 5;

    auto c1_test = bgv::EncryptionEngine::encrypt(s1_test, q_test, 1);
    auto c2_test = bgv::EncryptionEngine::encrypt(s1_test, q_test, 0);
    auto info_test = bgv::EncryptionEngine::encryptKeyPowersOf2(s2_test, q_test, s1_test);
    auto c_mult = bgv::EncryptionEngine::multiplyCyphertexts(c1_test, c1_test, info_test, q_test);
    std::cout << c1_test << "\n";
    std::cout << c2_test << "\n";
    std::cout << bgv::EncryptionEngine::decrypt(s1_test, q_test, c1_test) << "\n";
    std::cout << bgv::EncryptionEngine::decrypt(s1_test, q_test, c2_test) << "\n";
    std::cout <<  info_test << "\n";
    std::cout << c_mult << "\n";
    std::cout << bgv::EncryptionEngine::decrypt(s2_test, q_test, c_mult) << "\n";*/



    Eigen::VectorXi s(5);
    s << 3, 5, -2, 0, 1;

    Eigen::VectorXi s2(5);
    s2 << 7, 9, -1, 4, -15;
    int q = 65007;

    //std::cout << bgv::EncryptionEngine::createVectorBitDecomposition(s, 5) << "\n";
    //auto a = bgv::EncryptionEngine::encryptKeyPowersOf2(s2, q, s);
    //std::cout << a << "\n";

    for (int i = 0; i < 10 ;i ++)
    {
        auto c1 = bgv::EncryptionEngine::encrypt(s, q, 1);
        auto c2 = bgv::EncryptionEngine::encrypt(s, q, 0);
        auto c_sum = bgv::EncryptionEngine::addCyphertexts(c1, c2);
        auto m1 = bgv::EncryptionEngine::decrypt(s, q, c1);
        auto m2 = bgv::EncryptionEngine::decrypt(s, q, c2);
        auto m_sum = bgv::EncryptionEngine::decrypt(s, q, c_sum);

        auto info = bgv::EncryptionEngine::encryptKeyPowersOf2(s2, q, s);
        auto c_prod = bgv::EncryptionEngine::multiplyCyphertexts(c1, c1, info, q);
        auto m_prod = bgv::EncryptionEngine::decrypt(s2, q, c_prod);
        std::cout << m1 << "\n";
        std::cout << m2 << "\n";
        std::cout << m_sum << "\n";
        std::cout << m_prod << "\n";
        assert(m_prod(0, 0) == 1);
    }


    for (int i = 0; i < 10; i++)
    {
        auto c = bgv::EncryptionEngine::encrypt(s, q, 1);
        std::cout << c << "\n";
        auto m = bgv::EncryptionEngine::decrypt(s, q, c);
        std::cout << m << "\n";
        assert(m(0 , 0) == 1);
    }
    for (int i = 0; i < 10; i++)
    {
        auto c = bgv::EncryptionEngine::encrypt(s, q, 0);
        std::cout << c << "\n";
        auto m = bgv::EncryptionEngine::decrypt(s, q, c);
        std::cout << m << "\n";
        assert(m(0, 0) == 0);
    }
    return 0;
} 