#pragma once

#include "base/Exception.h"
#include <eigen3/Eigen/Core>

#ifdef HAVE_ARPACK
namespace NuTo
{
namespace ARPACKWRAP
{
extern "C" {
int dnaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol, double* resid, int* ncv, double* v,
            int* ldv, int* iparam, int* ipntr, double* workd, double* workl, int* lworkl, int* info);

int dneupd_(int* rvec, char* all, int* select, double* d_r, double* d_i, double* z, int* ldz, double* sigma_r,
            double* sigma_i, double* work_ev, char* bmat, int* n, char* which, int* nev, double* tol, double* resid,
            int* ncv, double* v, int* ldv, int* iparam, int* ipntr, double* workd, double* workl, int* lworkl,
            int* info);

int dsaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol, double* resid, int* ncv, double* v,
            int* ldv, int* iparam, int* ipntr, double* workd, double* workl, int* lworkl, int* info);

int dseupd_(int* rvec, char* all, int* select, double* d, double* z, int* ldz, double* sigma, char* bmat, int* n,
            char* which, int* nev, double* tol, double* resid, int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
            double* workd, double* workl, int* lworkl, int* info);
}
} // ARPACKWRAP
} // NuTo
#endif // HAVE_ARPACK


namespace NuTo
{
// forward declarations
template <class T>
class SparseMatrix;

namespace EIGEN_SOLVER_ARPACK
{
enum class eDriver;
enum class eWhich;
} // namespace EIGEN_SOLVER_ARPACK


//! @author Jörg F. Unger
//! @date Septemper 2013
//! @brief ... interface for the arpack eigenvalue solver
class EigenSolverArpack
{
public:
    //! @brief ... default constructor
    EigenSolverArpack();

    void Save (const std::string &filename, std::string rType )const
	{
		throw Exception("NuTo::EigenSolverArpack::Save] To be implemented.");
	}

	void Restore (const std::string &filename, std::string rType )
	{
		throw Exception("NuTo::EigenSolverArpack::Restore] To be implemented.");
	}

    //! @brief ... solve the eigenvalue problem of a single matrix
    void Solve(const NuTo::SparseMatrix<double>& rK, const NuTo::SparseMatrix<double>* rM, int rNumEigenValues,
               Eigen::MatrixXd& rEigenValues, Eigen::MatrixXd& rEigenVectors);

    //! @brief ... set the mode
    void SetDriver(NuTo::EIGEN_SOLVER_ARPACK::eDriver rDriver)
    {
        this->mDriver = rDriver;
    }

    //! @brief ... return mode
    NuTo::EIGEN_SOLVER_ARPACK::eDriver GetDriver() const
    {
        return this->mDriver;
    }

    //! @brief ... return mode
    NuTo::EIGEN_SOLVER_ARPACK::eWhich GetWhichEigenValues() const
    {
        return this->mWhich;
    }

    //! @brief ... set the which flag (LM Largest Magnitude, SM Smallest Magnitude, LA Largest Amplitude, SA Smallest
    //! Amplitude, BE
    void SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::eWhich rWhich)
    {
        this->mWhich = rWhich;
    }

    //! @brief ... set the mode
    void SetShift(double mSigmaR, double mSigmaI = 0.)
    {
        this->mSigmaR = mSigmaR;
        this->mSigmaI = mSigmaI;
    }

    //! @brief ... return shift real
    double GetShiftReal() const
    {
        return this->mSigmaR;
    }

    //! @brief ... return shift imag
    double GetShiftImag() const
    {
        return this->mSigmaI;
    }

    //! @brief calculates the largest eigenvalue + eigenvector using the driver 1 (normal mode)
    //! @param rM matrix
    //! @return pair<eValue, eVector>
    std::pair<double, Eigen::VectorXd> GetLargest(const SparseMatrix<double>& rM);

    //! @brief calculates the smallest eigenvalue + eigenvector using the driver 2 (shift inverted)
    //! @param rM matrix
    //! @return pair<eValue, eVector>
    std::pair<double, Eigen::VectorXd> GetSmallest(const SparseMatrix<double>& rM);

    bool GetShowTime()
    {
        return mShowTime;
    }

protected:
    //! @brief ... determines wether a refinement is enabled or disabled within the solution strategy
    NuTo::EIGEN_SOLVER_ARPACK::eDriver mDriver;
    NuTo::EIGEN_SOLVER_ARPACK::eWhich mWhich;
    double mTolerance; // tolerance, initially set to machine precision (relative eps)
    double mSigmaR; // real shift for Spectral transformation
    double mSigmaI; // imag shift for Spectral transformation
    bool mShowTime;
};
}
