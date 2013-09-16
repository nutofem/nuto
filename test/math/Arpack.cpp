#include <stdlib.h>
#include <boost/filesystem.hpp>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/math/EigenSolverArpack.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"

#define PRINTRESULT true

int main(int argc,char *argv[])
{
try
{
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> A_Full(8,8);
	A_Full(0,0) = 8.;
	A_Full(1,1) = 2.;
	A_Full(2,2) = 4.;
	A_Full(3,3) = 7.;
	A_Full(4,4) = 8.;
	A_Full(5,5) = 6.;
	A_Full(6,6) = 5.;
	A_Full(7,7) = 1.;
	A_Full(0,1) = 1.;
	A_Full(1,0) = 1.;
	A_Full(0,3) = -1.;
	A_Full(3,0) = -1.;
	A_Full(2,4) = 2.;
	A_Full(4,2) = 2.;
	A_Full(6,1) = -3.;
	A_Full(1,6) = -3.;

	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> M_Full(8,8);
	M_Full(0,0) = 1.;
	M_Full(1,1) = 1.;
	M_Full(2,2) = 1.;
	M_Full(3,3) = 1.;
	M_Full(4,4) = 1.;
	M_Full(5,5) = 1.;
	M_Full(6,6) = 1.;
	M_Full(7,7) = 1.;

	M_Full(0,1) = 0.5;
	M_Full(1,0) = 0.5;

	NuTo::SparseMatrixCSRVector2Symmetric<double> A_symmetric(A_Full);

	NuTo::SparseMatrixCSRVector2General<double> A_general(A_Full);

	NuTo::SparseMatrixCSRVector2Symmetric<double> M_symmetric(M_Full);

    if (PRINTRESULT)
	    std::cout << "FullMatrix\n" << A_Full << std::endl;
	//std::cout << "SymMatrix\n" << std::endl; A_symmetric.Info();
	//std::cout << "GenMatrix\n" <<  std::endl; A_general.Info();

    //values for the standard eigenvalue problem
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenVectorsStandard;
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValuesStandard;

    //values for the generalized eigenvalue problem
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenVectorsGeneral;
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValuesGeneral;

    int numEigenValuesCompute(3);

    //calculate eigenvalues of full matrix
    //std::cout << "calculate full matrix" << std::endl;
    A_Full.EigenVectorsSymmetric(eigenValuesStandard,eigenVectorsStandard);
    if (PRINTRESULT)
        std::cout << "eigenvalues of standard eigenvalue problem.\n" << eigenValuesStandard.Trans() << std::endl;
    A_Full.GeneralizedEigenVectorsSymmetric(M_Full,eigenValuesGeneral,eigenVectorsGeneral);
    if (PRINTRESULT)
        std::cout << "eigenvalues of generalized eigenvalue problem.\n" << eigenValuesGeneral.Trans() << std::endl;

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValuesStandardRef
    = eigenValuesStandard.GetBlock(8-numEigenValuesCompute,0,numEigenValuesCompute,1);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValuesGeneralRef
    = eigenValuesGeneral.GetBlock(8-numEigenValuesCompute,0,numEigenValuesCompute,1);

    //eigenvalue
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenVectors;
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValues;

    NuTo::EigenSolverArpack myEigenSolver;
    myEigenSolver.SetShowTime(false);


    //set driver for symmetric solution
    myEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::DSDRV1);
    myEigenSolver.SetShift(0.0); //not used
    //largest magnitude
    myEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::LA);
    myEigenSolver.Solve(A_symmetric,0,numEigenValuesCompute,eigenValues,eigenVectors);
    if (PRINTRESULT)
        std::cout << "largest eigenvalues of DSDRV1\n" << eigenValues.Trans() << std::endl;

    if ((eigenValues-eigenValuesStandardRef).norm()>1e-5)
    {
    	std::cout << "error calculating eigenvalues using DSDRV1 with error " << (eigenValues-eigenValuesStandardRef).norm() << std::endl;
        throw NuTo::MathException("Example Arpack : error in DSDRV1.");
        //std::cout << "eigenvectors\n" << eigenVectors << std::endl;
    }


    myEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::DSDRV2);
    myEigenSolver.SetShift(-10.0);
    //for DSDRV2, the which flag refers to the modified problem (A-shift*I)*x = lambda * x
    //the eigenvalues of the original problem are x_orig = shift + 1/x_mod
    //-->for a shift of zero, the largest eigenvalue is computed with SA (smallest amplitude)
    //if you want to calculate all eigenvalues closest to an given value, set Shift to this value and which to 'LM' (fabs - both directions)
    myEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::SA);
    myEigenSolver.Solve(A_symmetric,0,numEigenValuesCompute,eigenValues,eigenVectors);
    if (PRINTRESULT)
        std::cout << "largest eigenvalues of DSDRV2\n" << eigenValues.Trans() << std::endl;

    if ((eigenValues-eigenValuesStandardRef).norm()>1e-5)
    {
    	std::cout << "error calculating eigenvalues using DSDRV2 with error " << (eigenValues-eigenValuesStandardRef).norm() << std::endl;
        throw NuTo::MathException("Example Arpack : error in DSDRV2.");
        //std::cout << "eigenvectors\n" << eigenVectors << std::endl;
    }

    myEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::DSDRV3);
    myEigenSolver.SetShift(0.0); //not used
    //for DSDRV3, the which flag refers to the modified problem M^(-1)A*x = lambda * x
    myEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::LA);
    myEigenSolver.Solve(A_symmetric,&M_symmetric,numEigenValuesCompute,eigenValues,eigenVectors);
    if (PRINTRESULT)
        std::cout << "largest eigenvalues of DSDRV3\n" << eigenValues.Trans() << std::endl;

    if ((eigenValues-eigenValuesGeneralRef).norm()>1e-5)
    {
    	std::cout << "error calculating eigenvalues using DSDRV3 with error " << (eigenValues-eigenValuesGeneralRef).norm() << std::endl;
        throw NuTo::MathException("Example Arpack : error in DSDRV3.");
        //std::cout << "eigenvectors\n" << eigenVectors << std::endl;
    }

    myEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::DSDRV4);
    myEigenSolver.SetShift(-10.0);
    //for DSDRV4, the which flag refers to the modified problem (A-shift*M)^(-1) *M *x = lambda * x
    //the eigenvalues of the original problem are x_orig = shift + 1/x_mod
    //-->for a shift of zero, the largest eigenvalue is computed with SA (smallest amplitude)
    //if you want to calculate all eigenvalues closest to an given value, set Shift to this value and which to 'LM' (fabs - both directions)
    myEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::SA);
    myEigenSolver.Solve(A_symmetric,&M_symmetric,numEigenValuesCompute,eigenValues,eigenVectors);
    if (PRINTRESULT)
        std::cout << "largest eigenvalues of DSDRV4\n" << eigenValues.Trans() << std::endl;

    if ((eigenValues-eigenValuesGeneralRef).norm()>1e-5)
    {
    	std::cout << "error calculating eigenvalues using DSDRV4 with error " << (eigenValues-eigenValuesGeneralRef).norm() << std::endl;
        throw NuTo::MathException("Example Arpack : error in DSDRV4.");
        //std::cout << "eigenvectors\n" << eigenVectors << std::endl;
    }

    myEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::DSDRV5);
	myEigenSolver.SetShift(-10.0);
	myEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::LA);
	myEigenSolver.Solve(A_symmetric,&M_symmetric,numEigenValuesCompute,eigenValues,eigenVectors);
    if (PRINTRESULT)
        std::cout << "largest eigenvalues of DSDRV5\n" << eigenValues.Trans() << std::endl;

    if ((eigenValues-eigenValuesGeneralRef).norm()>1e-5)
    {
    	std::cout << "error calculating eigenvalues using DSDRV5 with error " << (eigenValues-eigenValuesGeneralRef).norm() << std::endl;
        throw NuTo::MathException("Example Arpack : error in DSDRV5.");
        //std::cout << "eigenvectors\n" << eigenVectors << std::endl;
    }

	myEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::DSDRV5);
	myEigenSolver.SetShift(-10.0);
	myEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::LA);
	myEigenSolver.Solve(A_symmetric,&M_symmetric,numEigenValuesCompute,eigenValues,eigenVectors);
    if (PRINTRESULT)
        std::cout << "largest eigenvalues of DSDRV6\n" << eigenValues.Trans() << std::endl;

    if ((eigenValues-eigenValuesGeneralRef).norm()>1e-5)
    {
    	std::cout << "error calculating eigenvalues using DSDRV6 with error " << (eigenValues-eigenValuesGeneralRef).norm() << std::endl;
        throw NuTo::MathException("Example Arpack : error in DSDRV6.");
        //std::cout << "eigenvectors\n" << eigenVectors << std::endl;
    }

    //set driver for unsymmetric solution
    myEigenSolver.SetDriver(NuTo::EIGEN_SOLVER_ARPACK::DNDRV1);
    myEigenSolver.SetShift(0.0); //not used
    //largest real
    myEigenSolver.SetWhichEigenValues(NuTo::EIGEN_SOLVER_ARPACK::LR);
    myEigenSolver.Solve(A_general,0,numEigenValuesCompute,eigenValues,eigenVectors);
    //sort the eigenvalues (I don't know why this is not done in ARPACK)
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValuesSorted = eigenValues.SortRow(0);
    if (PRINTRESULT)
        std::cout << "largest eigenvalues of DNDRV1\n" << eigenValuesSorted.Trans() << std::endl;

    if ((eigenValuesSorted.col(0)-eigenValuesStandardRef).norm()>1e-5)
    {
    	std::cout << "error calculating eigenvalues using DNDRV1 with error " << (eigenValuesSorted.col(0)-eigenValuesStandardRef).norm() << std::endl;
        throw NuTo::MathException("Example Arpack : error in DNDRV1.");
        //std::cout << "eigenvectors\n" << eigenVectors << std::endl;
    }

}
catch (NuTo::MathException& e)
{
    std::cout << "Error executing Arpack test routine "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    exit(-1);
}
catch (NuTo::Exception& e)
{
    std::cout << "Error executing Arpack "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    exit(-1);
}
catch (...)
{
    std::cout << "Error executing Arpack "<< std::endl;
    exit(-11);
}
}

