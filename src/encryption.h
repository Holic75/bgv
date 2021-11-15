#pragma once
#include "random_tools.h"
#include "zq_type.h"
#include <Eigen/Dense>

namespace bgv
{

/// \brief Set of utility functions for bgv encryption/decryption
/// 
/// TODO: optimize multiplication on "bit vectors" (i.e. the ones only containing -1, 1, 0)
/// TODO: account for short length of sk keyswitching info (if they are really short ?)
class EncryptionEngine
{

public:
    /// \brief Encrypt text with the specified key
	///
    /// @param s Secret key.
	/// @param modulus Modulus q of Zq.
	/// @param m text (vector of values of 0 or 1).
	/// @return Cyphertext corresponding to text m, every row(i) = encryption of m[i].
    static Eigen::MatrixXi encrypt(const Eigen::VectorXi& s, int modulus, const Eigen::VectorXi& m);

    /// \brief Encrypt a single bit value (0 or 1)
	///
    /// @param s Secret key.
	/// @param modulus Modulus q of Zq.
	/// @param m bit to encrypt (0 or 1).
	/// @return Encryption of bit.
    static Eigen::MatrixXi encrypt(const Eigen::VectorXi& s, int modulus, int value);


    /// \brief Compute side information for relinearization of old_s x old_s -> new_s
	///
    /// @param new_s New secret key.
	/// @param modulus Modulus q of Zq.
	/// @param old_s Old secret key.
	/// @return Information required for relinearization, every row contains encryption of 2^t s[i]s[j] indicies run as [t][j][i], stores only values for i>=j.
    static Eigen::MatrixXi encryptKeyPowersOf2(const Eigen::VectorXi& new_s, int modulus, const Eigen::VectorXi& old_s);

    /// \brief Compute side information for refresj of old_s x old_s -> new_s
	///
    /// @param new_s New secret key.
	/// @param new_modulus New modulus q of Zq.
	/// @param old_s Old secret key.
    /// @param old_modulus old modulus q of Zq.
	/// @return Information required for refresh, every row contains encryption of 2^k s[i]s[j]_t indicies run as [k][t][j][i], stores only values for i>=j.
    static Eigen::MatrixXi encryptKeyBitsPowersOf2(const Eigen::VectorXi& new_s, int new_modulus, const Eigen::VectorXi& old_s, int old_modulus);

    /// \brief Perform refresh
	///
    /// @param h Long expression: c = sum_ij h_ij s[i]s[j], stores only values for i >= j.
	/// @param info Info necessary for refresh (as provided by encryptKeyBitsPowersOf2)
    /// @param old_modulus Old modulus q of Zq.
    /// @param new_modulus New modulus q of Zq.
	/// @return Row vector: refreshed encoding c, corresponding to long expression with h.
    static Eigen::MatrixXi refresh(const Eigen::VectorXi& h, const Eigen::MatrixXi& info, int old_modulus, int new_modulus);

    /// \brief Compute long expression for multiplication of cyphertbits (?).
	///
    /// @param c1 Encryption of c1.
	/// @param c2 Encryption of c2 under the same key as c1.
    /// @param modulus Modulus q of Zq.
    /// @param new_modulus New modulus q of Zq.
	/// @return Long expression: c = sum_ij h_ij s[i]s[j], stores only values for i >= j.
    static Eigen::MatrixXi multiplyCyphertextsLong(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2, int modulus);

    /// \brief Compute product of cyphertbits (?).
	///
    /// @param c1 Encryption of c1.
	/// @param c2 Encryption of c2 under the same key as c1.
    /// @param info Info for relinearization as provided by encryptKeyPowersOf2.
    /// @param modulus Modulus q of Zq.
	/// @return encoding of product under new key from info.
    static Eigen::MatrixXi multiplyCyphertexts(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2, Eigen::MatrixXi info, int modulus);


    /// \brief Relinearize long expression 
	///
    /// @param h_vec  Elements of h_ij, where c = sum_ij h_ij s[i]s[j], stores only values for i >= j.
    /// @param info Info for relinearization as provided by encryptKeyPowersOf2.
	/// @return encoding of long expression under new key from info.
    static Eigen::MatrixXi relinearize(const Eigen::VectorXi& h_vec,  Eigen::MatrixXi info);

    /// \brief Directly compute sum of cyphertbits (?).
	///
    /// @param c1 Encryption of c1.
	/// @param c2 Encryption of c2 under the same key as c1.
    /// @param modulus Modulus q of Zq.
	/// @return encoding of sum under the same key as c1 and c2.
    static Eigen::MatrixXi addCyphertexts(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2, int modulus);

    /// \brief Obtain long expression corresponding to the sum of c1 and c2.
	///
    /// @param c1 Encryption of c1.
	/// @param c2 Encryption of c2 under the same key as c1.
    /// @param modulus Modulus q of Zq.
	/// @return Long expression: c = sum_ij h_ij s[i]s[j], stores only values for i >= j.
    static Eigen::VectorXi addCyphertextsLong(const Eigen::MatrixXi& c1, const Eigen::MatrixXi& c2, int modulus);

    /// \brief Decrypt cyphertext
    /// @param s Secret key.
	/// @param modulus Modulus q of Zq.
	/// @param cyphertext Cyphertext encoded under key s.
	/// @return Text corresponding to cyphertext
    static Eigen::VectorXi decrypt(const Eigen::VectorXi& s, int modulus, const Eigen::MatrixXi& cyphertext);

    static int log2(unsigned int q)
    {
        int count = 0;
        while (q > 0)
        {
            count++;
            q = q >> 1;
        }
        return count;
    }

private:
    static Eigen::VectorXi storeSymmetricMatrixAsVector(const Eigen::MatrixXi& m, bool sum);
    static Eigen::VectorXi createVectorBitDecomposition(const Eigen::VectorXi& v, int n_bits);
    static Eigen::VectorXi createPowersOf2Vector(const Eigen::MatrixXi& vec, int n_bits);
    static  Eigen::VectorXi createPowersOf2(int size);
};

}