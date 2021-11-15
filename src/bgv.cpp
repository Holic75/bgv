#include "zq_type.h"
#include "random_tools.h"
#include "encryption.h"
#include <Eigen/Dense>
#include <vector>
#include <utility>

int main (int argc, char *argv[]) 
{
    //keys, normally should be generated randomly
    //but we will use fixed values here for simplicity
    Eigen::VectorXi s(5);
    s << 3, 5, -2, 0, 1;

    Eigen::VectorXi s2(5);
    s2 << 7, 9, -1, 4, -15;

    int q = 32771;
    int q2 = 16273;


    std::vector<std::pair<int, int>> bits(4);
    bits[0] = {0, 0};
    bits[1] = {1, 0};
    bits[2] = {0, 1};
    bits[3] = {1, 1};

    for (auto& b: bits)
    {
        std::cout << "Starting test\n";
        std::cout << b.first << " <-bit1\n";
        std::cout << b.second << " <-bit2\n";

        auto c1 = bgv::EncryptionEngine::encrypt(s, q, b.first);
        auto c2 = bgv::EncryptionEngine::encrypt(s, q, b.second);       
        auto m1 = bgv::EncryptionEngine::decrypt(s, q, c1);
        auto m2 = bgv::EncryptionEngine::decrypt(s, q, c2);

        std::cout << "Test validity of encryption / decryption\n";
        std::cout << m1 << " <-m1\n";
        std::cout << m2 << " <-m2\n";
        assert(m1(0, 0) == b.first);
        assert(m2(0, 0) == b.second);

        auto c_sum = bgv::EncryptionEngine::addCyphertexts(c1, c2, q);
        auto m_sum = bgv::EncryptionEngine::decrypt(s, q, c_sum);

        std::cout << "Test validity of addition\n";
        int expected_sum = (b.first + b.second) % 2;
        std::cout << expected_sum  << "<- bit1 + bit2\n";
        std::cout << m_sum(0, 0) << "<- m1 + m2\n";
        assert(m_sum(0, 0) == expected_sum);

        auto info = bgv::EncryptionEngine::encryptKeyPowersOf2(s2, q, s);
        auto c_prod = bgv::EncryptionEngine::multiplyCyphertexts(c1, c2, info, q);
        auto m_prod = bgv::EncryptionEngine::decrypt(s2, q, c_prod);

        std::cout << "Test validity of multiplication\n";
        int expected_prod = b.first * b.second;
        std::cout << expected_prod  << "<- bit1 * bit2\n";
        std::cout << m_prod(0, 0) << "<- m1 * m2\n";
        assert(m_prod(0, 0) == expected_prod);

        auto c_sum_long = bgv::EncryptionEngine::addCyphertextsLong(c1, c2, q);
        auto c_sum_relin = bgv::EncryptionEngine::relinearize(c_sum_long, info);
        auto m_sum_relin = bgv::EncryptionEngine::decrypt(s2, q, c_sum_relin);

        std::cout << "Test addition with relinearization\n";
        std::cout << expected_sum  << "<- bit1 + bit2\n";
        std::cout << m_sum_relin(0, 0) << "<- m1 + m2 (relinearization)\n";
        assert(m_sum_relin(0, 0) == expected_sum);


        auto info_sm = bgv::EncryptionEngine::encryptKeyBitsPowersOf2(s2, q2, s, q);
        auto c_sum_refresh = bgv::EncryptionEngine::refresh(c_sum_long, info_sm, q, q2);
        auto m_sum_refresh = bgv::EncryptionEngine::decrypt(s2, q2, c_sum_refresh);

        std::cout << "Test addition with refresh\n";
        std::cout << expected_sum  << "<- bit1 + bit2\n";
        std::cout << m_sum_refresh(0, 0) << "<- m1 + m2 (refresh)\n";
        assert(m_sum_refresh(0, 0) == expected_sum);

        auto c_prod_long = bgv::EncryptionEngine::multiplyCyphertextsLong(c1, c2, q);
        auto c_prod_refresh = bgv::EncryptionEngine::refresh(c_prod_long, info_sm, q, q2);
        auto m_prod_refresh = bgv::EncryptionEngine::decrypt(s2, q2, c_prod_refresh);

        std::cout << "Test multiplication with refresh\n";
        std::cout << expected_prod  << "<- bit1 * bit2\n";
        std::cout << m_prod_refresh(0, 0) << "<- m1 * m2 (refresh)\n";
        assert(m_prod_refresh(0, 0) == expected_prod);

        std::cout << "Test ok\n\n";
    }

    return 0;
} 