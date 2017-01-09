#include <limits>       // std::numeric_limits
#include "math/EigenSolverArpack.h"
#include "math/EigenSolverArpackEnum.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/MathException.h"
#include "base/Timer.h"
#include <memory>

NuTo::EigenSolverArpack::EigenSolverArpack() : NuToObject()
{
#ifdef HAVE_ARPACK
	mDriver  = NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV1;
	mWhich = NuTo::EIGEN_SOLVER_ARPACK::eWhich::LM;
    mTolerance = std::numeric_limits<double>::epsilon();//machine precision
    mSigmaR = 0.; //real shift used for spectral transformations
    mSigmaI = 0.; //imag shift used for spectral transformations
#else
    throw MathException(__PRETTY_FUNCTION__, "NuTo wasn't compiled with ARPACK.");
#endif
}

void NuTo::EigenSolverArpack::Solve(const NuTo::SparseMatrix<double>& rK,
		const NuTo::SparseMatrix<double>* rM,
		int rNumEigenValues,
        Eigen::MatrixXd& rEigenValues,
        Eigen::MatrixXd& rEigenVectors)
{
#ifdef HAVE_ARPACK
    NuTo::Timer(__FUNCTION__, GetShowTime());
    char whichEigenValue[3];  //LA SA SM BE

    if (rK.IsSymmetric())
    {
	    //use dsaupd_ and dseupd_
    	switch (mWhich)
		{
        case NuTo::EIGEN_SOLVER_ARPACK::eWhich::LA:
			whichEigenValue[0] = 'L';
			whichEigenValue[1] = 'A';
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eWhich::SA:
			whichEigenValue[0] = 'S';
			whichEigenValue[1] = 'A';
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eWhich::LM:
			whichEigenValue[0] = 'L';
			whichEigenValue[1] = 'M';
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eWhich::SM:
			whichEigenValue[0] = 'S';
			whichEigenValue[1] = 'M';
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eWhich::BE:
			whichEigenValue[0] = 'B';
			whichEigenValue[1] = 'E';
			break;
		default:
			throw MathException(__PRETTY_FUNCTION__, "which type not implemented for symmetric matrices (LA, SA, LM, SM, BE.");
		}
    }
    else
    {
	    //use dnaupd_ and dneupd_
    	switch (mWhich)
		{
		case NuTo::EIGEN_SOLVER_ARPACK::eWhich::LM:
			whichEigenValue[0] = 'L';
			whichEigenValue[1] = 'M';
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eWhich::SM:
			whichEigenValue[0] = 'S';
			whichEigenValue[1] = 'M';
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eWhich::LR:
			whichEigenValue[0] = 'L';
			whichEigenValue[1] = 'R';
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eWhich::SR:
			whichEigenValue[0] = 'S';
			whichEigenValue[1] = 'R';
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eWhich::LI:
			whichEigenValue[0] = 'L';
			whichEigenValue[1] = 'I';
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eWhich::SI:
			whichEigenValue[0] = 'S';
			whichEigenValue[1] = 'I';
			break;
		default:
			throw MathException(__PRETTY_FUNCTION__, "which type not implemented for general matrices (LM, SM, LR, SR, LI, SI).");
		}
    }
    if (GetVerboseLevel()>5)
        std::cout << "machine tolerance " << mTolerance << std::endl;

    int ido=0; //return value for the reverse communication loop
    char bmat; //standard 'I' or general 'G' eigenvalue problem
    int n = rK.GetNumRows();
    if (n!=rK.GetNumColumns())
    	throw MathException(__PRETTY_FUNCTION__, "only square matrices are supported");
    //mWhichEigenValue[2];  //LA SA SM BE w
    int nev = rNumEigenValues;
    if (rNumEigenValues>=n-1)
    	throw MathException(__PRETTY_FUNCTION__, "the number of extracted eigenvalues must be less smaller than n-1 (n:dimension of matrix)");

    double tol_err = mTolerance; //solution tolerance (stopping criterion)
    std::vector<double> resid(n);  //residual vector
    int ncv=std::min(2*nev+1, n);    //number of Arnoldi vectors used, describes the size of array v
    //std::cout << "ncv " << ncv << " nev " << nev << std::endl;
    std::vector<double> v(n*ncv);
    int ldv=n; //leading dimension of v

    int iParam[11];
    iParam[0]  = 1; //1: exact shifts  0: given by the user via reverse communication
    iParam[1]  = 0; //not used
    iParam[2]  = 500; //input: maximum number of Arnoldi update iterations, output: actual number used (dimension of Krylow subspace)
    iParam[3]  = 1; //blocksize, has to be set to 1
    iParam[4]  = 0; //number of converged Ritz values (eigenvalues) that fulfill the convergence property)
    iParam[5]  = 0; //not used
    //iParam[6]  = 1; // set in a following switch case mDriver
    iParam[7]  = 0; //used as output to the user for iParam==0, gives the number of shifts the user has to provide
    iParam[8]  = 0; //output: total number of OP*x operations (matrix * vector)
    iParam[9]  = 0; //output: total number of B*x operations (matrix * vector)
    iParam[10] = 0; //output: total number of reorthogonalization steps

    int iPntr[14]; //output : Pointer to mark the starting locations in the WORKD and WORKL
                   //arrays for matrices/vectors used by the Arnoldi iteration.

    FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> workd(n,3);
    int lworkl;

    int info(0);

    bool solverRequired(false);
    std::unique_ptr<NuTo::SparseMatrixCSR<double>> solveMatrix = nullptr;
    Eigen::VectorXd solution;

    SparseDirectSolverMUMPS solver;
    solver.SetShowTime(GetShowTime());

    if (rK.IsSymmetric())
    {
    	lworkl=ncv*(ncv+8);
		switch (this->mDriver)
		{
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV1:
			bmat = 'I';
			iParam[6] = 1; //dsdrv1
			solverRequired = false;
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV2:
			bmat = 'I';
			iParam[6] = 3; //dsdrv2
			solverRequired = true;
			switch (rK.GetSparseMatrixType())
			{
			case eSparseMatrixType::CSRGENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRGENERAL to SparseMatrixCSRSymmetric not possible.");
			break;
			case eSparseMatrixType::CSRSYMMETRIC:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRSymmetric<double>>(rK.AsSparseMatrixCSRSymmetric());
			break;
			case eSparseMatrixType::CSRVECTOR2GENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRVECTOR2GENERAL to SparseMatrixCSRSymmetric not possible.");
			break;
			case eSparseMatrixType::CSRVECTOR2SYMMETRIC:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRSymmetric<double>>(rK.AsSparseMatrixCSRVector2Symmetric());
			break;
			default:
		    	throw MathException(__PRETTY_FUNCTION__, "matrix type not implemented.");
			};
			//add diagonal shift
            for (int count=0; count<solveMatrix->GetNumRows(); count++)
			{
				solveMatrix->AddValue(count,count,-mSigmaR);
			}
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV3:
			bmat = 'G';
			iParam[6] = 2; //dsdrv3
			solverRequired = true;
			if (rM==0)
				throw MathException(__PRETTY_FUNCTION__, "second matrix (M) required.");
			switch (rM->GetSparseMatrixType())
			{
			case eSparseMatrixType::CSRGENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRGENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRSYMMETRIC:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRSymmetric<double>>(rM->AsSparseMatrixCSRSymmetric());
			break;
			case eSparseMatrixType::CSRVECTOR2GENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRVECTOR2GENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRVECTOR2SYMMETRIC:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRSymmetric<double>>(rM->AsSparseMatrixCSRVector2Symmetric());
			break;
			default:
		    	throw MathException(__PRETTY_FUNCTION__, "matrix type not implemented.");
			};
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV4:
			bmat = 'G';
			iParam[6] = 3; //dsdrv4
			solverRequired = true;
			if (rM==0)
				throw MathException(__PRETTY_FUNCTION__, "second matrix (M) required.");

			switch (rK.GetSparseMatrixType())
			{
			case eSparseMatrixType::CSRGENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRGENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRSYMMETRIC:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRSymmetric<double>>(rK.AsSparseMatrixCSRSymmetric());
			break;
			case eSparseMatrixType::CSRVECTOR2GENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRVECTOR2GENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRVECTOR2SYMMETRIC:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRSymmetric<double>>(rK.AsSparseMatrixCSRVector2Symmetric());
			break;
			default:
		    	throw MathException(__PRETTY_FUNCTION__, "matrix type not implemented.");
			};
			//subtract shift
			switch (rM->GetSparseMatrixType())
			{
			case eSparseMatrixType::CSRGENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRGENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRSYMMETRIC:
				(*solveMatrix) += rM->AsSparseMatrixCSRSymmetric()*(-mSigmaR);
			break;
			case eSparseMatrixType::CSRVECTOR2GENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRVECTOR2GENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRVECTOR2SYMMETRIC:
				(*solveMatrix) += rM->AsSparseMatrixCSRVector2Symmetric()*(-mSigmaR);
			break;
			default:
		    	throw MathException(__PRETTY_FUNCTION__, "matrix type not implemented.");
			};
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV5:
			bmat = 'G';
			iParam[6] = 4; //dsdrv5
			solverRequired = true;
			if (rM==0)
				throw MathException(__PRETTY_FUNCTION__, "second matrix (M) required.");

			switch (rK.GetSparseMatrixType())
			{
			case eSparseMatrixType::CSRGENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRGENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRSYMMETRIC:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRSymmetric<double>>(rK.AsSparseMatrixCSRSymmetric());
			break;
			case eSparseMatrixType::CSRVECTOR2GENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRVECTOR2GENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRVECTOR2SYMMETRIC:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRSymmetric<double>>(rK.AsSparseMatrixCSRVector2Symmetric());
			break;
			default:
		    	throw MathException(__PRETTY_FUNCTION__, "matrix type not implemented.");
			};
			//subtract shift
			switch (rM->GetSparseMatrixType())
			{
			case eSparseMatrixType::CSRGENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRGENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRSYMMETRIC:
				(*solveMatrix) += rM->AsSparseMatrixCSRSymmetric()*(-mSigmaR);
			break;
			case eSparseMatrixType::CSRVECTOR2GENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRVECTOR2GENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRVECTOR2SYMMETRIC:
				(*solveMatrix) += rM->AsSparseMatrixCSRVector2Symmetric()*(-mSigmaR);
			break;
			default:
		    	throw MathException(__PRETTY_FUNCTION__, "matrix type not implemented.");
			};
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV6:
			bmat = 'G';
			iParam[6] = 5; //dsdrv6
			solverRequired = true;
			if (rM==0)
				throw MathException(__PRETTY_FUNCTION__, "second matrix (M) required.");

			switch (rK.GetSparseMatrixType())
			{
			case eSparseMatrixType::CSRGENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRGENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRSYMMETRIC:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRSymmetric<double>>(rK.AsSparseMatrixCSRSymmetric());
			break;
			case eSparseMatrixType::CSRVECTOR2GENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRVECTOR2GENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRVECTOR2SYMMETRIC:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRSymmetric<double>>(rK.AsSparseMatrixCSRVector2Symmetric());
			break;
			default:
		    	throw MathException(__PRETTY_FUNCTION__, "matrix type not implemented.");
			};
			//subtract shift
			switch (rM->GetSparseMatrixType())
			{
			case eSparseMatrixType::CSRGENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRGENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRSYMMETRIC:
				(*solveMatrix) += rM->AsSparseMatrixCSRSymmetric()*(-mSigmaR);
			break;
			case eSparseMatrixType::CSRVECTOR2GENERAL:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRVECTOR2GENERAL to SparseMatrixCSRSymmetric not implemented.");
			break;
			case eSparseMatrixType::CSRVECTOR2SYMMETRIC:
				(*solveMatrix) += rM->AsSparseMatrixCSRVector2Symmetric()*(-mSigmaR);
			break;
			default:
		    	throw MathException(__PRETTY_FUNCTION__, "matrix type not implemented.");
			};
			break;
		default:
			throw MathException(__PRETTY_FUNCTION__, "use a symmetric driver DSDRV1 .. DSDRV6");
		}
    }
    else
    {
    	lworkl=3*ncv*(ncv+6);
		switch (this->mDriver)
		{
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DNDRV1:
			bmat = 'I';
			iParam[6] = 1; //dndrv1
			solverRequired = false;
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DNDRV2:
            bmat = 'I';
			iParam[6] = 3; //dndrv2
			solverRequired = true;
			switch (rK.GetSparseMatrixType())
			{
			case eSparseMatrixType::CSRGENERAL:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRGeneral<double>>(rK.AsSparseMatrixCSRGeneral());
			break;
			case eSparseMatrixType::CSRSYMMETRIC:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRSYMMETRIC to SparseMatrixCSRGeneral not possible.");
			break;
			case eSparseMatrixType::CSRVECTOR2GENERAL:
				solveMatrix = std::make_unique<NuTo::SparseMatrixCSRGeneral<double>>(rK.AsSparseMatrixCSRVector2General());
			break;
			case eSparseMatrixType::CSRVECTOR2SYMMETRIC:
				throw MathException(__PRETTY_FUNCTION__, "Conversion from CSRVECTOR2SYMMETRIC to SparseMatrixCSRGeneral not possible.");
			break;
			default:
		    	throw MathException(__PRETTY_FUNCTION__, "matrix type not implemented.");
			};
			//add diagonal shift
            for (int count=0; count<solveMatrix->GetNumRows(); count++)
			{
				solveMatrix->AddValue(count,count,-mSigmaR);
			}
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DNDRV3:
			bmat = 'G';
			iParam[6] = 2; //dndrv3
			solverRequired = true;
			throw MathException(__PRETTY_FUNCTION__, "Driver DNDRV3 not implemented.");
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DNDRV4:
			bmat = 'G';
			iParam[6] = 4; //dndrv4
			solverRequired = true;
			throw MathException(__PRETTY_FUNCTION__, "Driver DNDRV4 not implemented.");
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DNDRV5:
			bmat = 'G';
			iParam[6] = 3; //dndrv5
			solverRequired = true;
			throw MathException(__PRETTY_FUNCTION__, "Driver DNDRV5 not implemented.");
			break;
		case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DNDRV6:
			bmat = 'G';
			iParam[6] = 4; //dndrv6
			solverRequired = true;
			throw MathException(__PRETTY_FUNCTION__, "Driver DNDRV6 not implemented.");
			break;
		default:
			throw MathException(__PRETTY_FUNCTION__, "use an unsymmetric driver DNDRV1 .. DNDRV6");
		}
    }
    std::vector<double> workl(lworkl);

    //initialize solver (factorization of M) for all other than the regular mode
    if (solverRequired)
    {
    	solveMatrix->SetOneBasedIndexing();
    	solver.Factorization(*solveMatrix);
        solver.SetShowTime(false);
    }

    do
    {
    	if (rK.IsSymmetric())
    	{
            ARPACKWRAP::dsaupd_(&ido, &bmat, &n, whichEigenValue, &nev,
                                &tol_err, &resid[0], &ncv, &v[0], &ldv,
                                iParam, iPntr, workd.data(), &workl[0],
                                &lworkl, &info);
    	}
    	else
    	{
            ARPACKWRAP::dnaupd_ (&ido, &bmat, &n, whichEigenValue, &nev,
                                 &tol_err, &resid[0], &ncv, &v[0], &ldv,
                                 iParam, iPntr,  workd.data(), &workl[0],
                                 &lworkl, &info);
    	}

    	if (info!=0)
    	{
    		switch (info)
    		{
    		case 1:
    			std::cout << "number of converged Ritz values " << iParam[4] << std::endl;
    			std::cout << "Ritz values\n real:" << Eigen::Map<Eigen::VectorXd>(&(workl[iPntr[5]-1]),iParam[4]).transpose()<< std::endl;
    			throw MathException(__PRETTY_FUNCTION__, "Error calling dnaupd_ (1) maximum number of iterations performed.");

    			break;
    		case -3:
    			throw MathException(__PRETTY_FUNCTION__, "Error calling dnaupd_ (-3) NCV must be greater than NEW and less then or equal to N.");
    			break;
    		case -4:
    			throw MathException(__PRETTY_FUNCTION__, "Error calling dnaupd_ (-4) the maximum number of Arnoldi iterations must be positive.");
    			break;
    		case -5:
    			throw MathException(__PRETTY_FUNCTION__, "Error calling dnaupd_ (-5) which must be eitehr LM SM, LA SA or BE.");
    			break;

    		default:
    			std::cout << "info " << info << std::endl;
    			throw MathException(__PRETTY_FUNCTION__, "Error calling dnaupd_");
    		}
    	}

        switch (this->mDriver)
        {
        case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV1:
        case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DNDRV1:
        	switch(ido)
        	{
        	case 1:
        	case -1:
        		//matrix multiplication
        		Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = rK.operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n));
        	break;
        	case 99:
       		break;
        	default:
            	throw MathException(__PRETTY_FUNCTION__, "ido flag of dnaupd_/dsaupd_ not supported.");
        	}
        	break;
        case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV2:
        case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DNDRV2:
        	switch(ido)
        	{
        	case 1:
        	case -1:
        		//solve OP*w = v
        		solver.Solution(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n), solution);
        		Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = solution;
        	break;
        	case 99:
       		break;
        	default:
            	throw MathException(__PRETTY_FUNCTION__, "ido flag of dnaupd_/dsaupd_ not supported.");
        	}
        	break;
        case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV3:
        	switch(ido)
        	{
        	case 1:
        	case -1:
        		//compute A*v = y (store temporary in pos 1
                Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = rK.operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n));
                //store y = pos 0 (required)
                Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n) = Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n);
        		//solve M*w = y
        		solver.Solution(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n), solution);
        		Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = solution;
        	break;
        	case 2:
        		//compute M*v = w
                Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = rM->operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n));
        	break;
        	case 99:
       		break;
        	default:
            	throw MathException(__PRETTY_FUNCTION__, "ido flag of dnaupd_/dsaupd_ not supported.");
        	}
        	break;
        case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV4:
        	switch(ido)
        	{
        	case -1:
        		//compute M*v = y (store temporary in pos 1
                Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = rM->operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n));
        		//solve OP * w = y
        		solver.Solution(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n), solution);
        		Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = solution;
            break;
        	case 1:
        		//solve OP * w = y
        		solver.Solution(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[2]-1]),n), solution);
        		Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = solution;
        	break;
        	case 2:
        		//compute M*v = w
                Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = rM->operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n));
        	break;
        	case 99:
       		break;
        	default:
            	throw MathException(__PRETTY_FUNCTION__, "ido flag of dnaupd_/dsaupd_ not supported.");
        	}
        	break;
        case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV5:
        	switch(ido)
        	{
        	case -1:
        		//compute M*v = y (store temporary in pos 1
                Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = rK.operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n));
        		//solve OP * w = y
        		solver.Solution(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n), solution);
        		Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = solution;
            break;
        	case 1:
        		//solve OP * w = y
        		solver.Solution(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[2]-1]),n), solution);
        		Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = solution;
        	break;
        	case 2:
        		//compute M*v = w
                Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = rK.operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n));
        	break;
        	case 99:
       		break;
        	default:
            	throw MathException(__PRETTY_FUNCTION__, "ido flag of dnaupd_/dsaupd_ not supported.");
        	}
        	break;
        case NuTo::EIGEN_SOLVER_ARPACK::eDriver::DSDRV6:
        	switch(ido)
        	{
        	case -1:
        		//compute (A+mSigmaR*M)*v = y (store temporary in pos 1
                Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = rK.operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n))
                + rM->operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n)*mSigmaR);
        		//solve OP * w = y
        		solver.Solution(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n), solution);
        		Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = solution;
            break;
        	case 1:
        		//compute (A+mSigmaR*M)*v = y (store temporary in pos 1
                Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = rK.operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[2]-1]),n))
                + rM->operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n)*mSigmaR);
        		//solve OP * w = y
        		solver.Solution(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n), solution);
        		Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = solution;
        	break;
        	case 2:
        		//compute M*v = w
                Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[1]-1]),n) = rM->operator*(Eigen::Map<Eigen::VectorXd>(&(workd.data()[iPntr[0]-1]),n));
        	break;
        	case 99:
       		break;
        	default:
            	throw MathException(__PRETTY_FUNCTION__, "ido flag of dnaupd_/dsaupd_ not supported.");
        	}
        	break;
        default:
        	throw MathException(__PRETTY_FUNCTION__, "mode/driver not correctly implemented.");
        }
    }
    while (ido!=99);

    if (info==0)
    {
    	std::vector<int> select(ncv);
    	int rVec(1); //0: don't calculate eigenvectors or 1: calculate eigenvectors
    	char howmany('A');
    	if (rK.IsSymmetric())
    	{
			rEigenValues.Resize(rNumEigenValues,1);
			rEigenVectors.Resize(n,rNumEigenValues);
			int ldz(n);
            ARPACKWRAP::dseupd_ (&rVec, &howmany, &select[0], rEigenValues.col(0).data(),
                                 rEigenVectors.data(), &ldz,	&mSigmaR,
                                 &bmat, &n, whichEigenValue, &nev,
                                 &tol_err, &resid[0], &ncv, &v[0], &ldv,
                                 iParam, iPntr,  workd.data(), &workl[0],
                                 &lworkl, &info);
	    	if (info!=0)
	    		throw MathException(__PRETTY_FUNCTION__, "Error calling dseupd_.");
   	    }
    	else
    	{
			rEigenValues.Resize(rNumEigenValues,2);
			rEigenVectors.Resize(n,rNumEigenValues);
			int ldz(n);
			std::vector<double> workev(3*ncv);
            ARPACKWRAP::dneupd_ (&rVec, &howmany, &select[0], rEigenValues.col(0).data(), rEigenValues.col(1).data(),
                                 rEigenVectors.data(), &ldz,
                                 &mSigmaR, &mSigmaI, &workev[0],
                                 &bmat, &n, whichEigenValue, &nev,
                                 &tol_err, &resid[0], &ncv, &v[0], &ldv,
                                 iParam, iPntr,  workd.data(), &workl[0],
                                 &lworkl, &info);
	    	if (info!=0)
	    		throw MathException(__PRETTY_FUNCTION__, "Error calling dneupd_.");
    	}

    }
    else
    {
		throw MathException(__PRETTY_FUNCTION__, "Error calculating eigenvalues.");
    }
#else //HAVE_ARPACK
    throw MathException(__PRETTY_FUNCTION__, "NuTo wasn't compiled with ARPACK.");
#endif//HAVE_ARPACK
}


std::pair<double, Eigen::VectorXd> NuTo::EigenSolverArpack::GetLargest(const SparseMatrix<double>& rM)
{
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValues;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenVectors;

    if (rM.IsSymmetric())
    {
        mDriver = EIGEN_SOLVER_ARPACK::eDriver::DSDRV1;
        mWhich = EIGEN_SOLVER_ARPACK::eWhich::LM;
    }
    else
    {
        mDriver = EIGEN_SOLVER_ARPACK::eDriver::DNDRV1;
        mWhich = EIGEN_SOLVER_ARPACK::eWhich::LR;
    }
    Solve(rM, nullptr, 1, eigenValues, eigenVectors);
    return std::make_pair(eigenValues(0,0), eigenVectors.col(0));
}

std::pair<double, Eigen::VectorXd> NuTo::EigenSolverArpack::GetSmallest(const SparseMatrix<double>& rM)
{
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenValues;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenVectors;

    if (rM.IsSymmetric())
    {
        mDriver = EIGEN_SOLVER_ARPACK::eDriver::DSDRV2;
        mWhich = EIGEN_SOLVER_ARPACK::eWhich::LM;
    }
    else
    {
        mDriver = EIGEN_SOLVER_ARPACK::eDriver::DNDRV2;
        mWhich = EIGEN_SOLVER_ARPACK::eWhich::LR;
    }
    Solve(rM, nullptr, 2, eigenValues, eigenVectors);
    return std::make_pair(eigenValues(0,0), eigenVectors.col(0));
}
