// $Id: SparseDirectSolverMKLDSS.h 625 2013-04-22 16:37:11Z unger3 $

#pragma once


#include "nuto/base/NuToObject.h"
#include "nuto/math/MathException.h"
#include <eigen3/Eigen/Core>


namespace NuTo
{
// forward declarations
template<class T> class SparseMatrix;
template<class T, int rows> class FullVector;
template <class T, int rows, int cols> class FullMatrix;

namespace EIGEN_SOLVER_ARPACK
{
    enum class eDriver;
    enum class eWhich;
}// namespace EIGEN_SOLVER_ARPACK


//! @author Jörg F. Unger
//! @date Septemper 2013
//! @brief ... interface for the arpack eigenvalue solver
class EigenSolverArpack  : public NuToObject
{
public:
    //! @brief ... default constructor
    EigenSolverArpack();

    //! @brief ... print information about the class attributes
    void Info()const
    {
    }

    void Save (const std::string &filename, std::string rType )const
	{
		throw MathException("NuTo::EigenSolverArpack::Save] To be implemented.");
	}

	void Restore (const std::string &filename, std::string rType )
	{
		throw MathException("NuTo::EigenSolverArpack::Restore] To be implemented.");
	}

    std::string GetTypeId()const
    {
        return std::string("EigenSolverArpack");
    }

    //! @brief ... solve the eigenvalue problem of a single matrix
    void Solve(const NuTo::SparseMatrix<double>& rK,
    		const NuTo::SparseMatrix<double>* rM,
    		int rNumEigenValues,
    		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEigenValues,
    		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEigenVectors);

    //! @brief ... set the mode
    inline void SetDriver(NuTo::EIGEN_SOLVER_ARPACK::eDriver rDriver)
    {
        this->mDriver = rDriver;
    }

    //! @brief ... return mode
    inline NuTo::EIGEN_SOLVER_ARPACK::eDriver GetDriver()const
    {
        return this->mDriver;
    }

    //! @brief ... return mode
    inline NuTo::EIGEN_SOLVER_ARPACK::eWhich GetWhichEigenValues()const
    {
        return this->mWhich;
    }

    //! @brief ... set the which flag (LM Largest Magnitude, SM Smallest Magnitude, LA Largest Amplitude, SA Smallest Amplitude, BE
    inline void SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::eWhich rWhich)
    {
        this->mWhich = rWhich;
    }

    //! @brief ... set the mode
    inline void SetShift(double mSigmaR, double mSigmaI=0.)
    {
        this->mSigmaR = mSigmaR;
        this->mSigmaI = mSigmaI;
    }

    //! @brief ... return shift real
    inline double GetShiftReal()const
    {
        return this->mSigmaR;
    }

    //! @brief ... return shift imag
    inline double GetShiftImag()const
    {
        return this->mSigmaI;
    }

protected:
    //! @brief ... determines wether a refinement is enabled or disabled within the solution strategy
    NuTo::EIGEN_SOLVER_ARPACK::eDriver mDriver;
    NuTo::EIGEN_SOLVER_ARPACK::eWhich mWhich;
    double mTolerance; //tolerance, initially set to machine precision (relative eps)
    double mSigmaR; //real shift for Spectral transformation
    double mSigmaI; //imag shift for Spectral transformation

};
}
