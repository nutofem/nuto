//! @author Andrea Keßler, ISM
//! @brief ... von-Mises-Wielandt method to calculate max and min eigenvalues and vectors

#include <iostream>
#include <fstream>
#include "optimize/MisesWielandt.h"
#include "base/Timer.h"
#define machine_precision 1e-15
//sqrt machine_precision

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)

int NuTo::MisesWielandt::Optimize()
{
    NuTo::Timer timer(__PRETTY_FUNCTION__, mShowTime);
    
    const double tol=mAccuracyGradient;
	double lambda=0.,		//lambda_min = lambda
		   prevLambda=0.,
           norm=0.;
           //prevNorm=0.;

    int numGradientCalls(0);   // number of gradient calls

    eOptimizationReturnAttributes returnValue;

	std::fstream outputTime;
	std::string filename = "timeOutput";

	// u  = parameters
	// r  = gradient

	std::vector<double> u(mNumParameters+3,1.);
	std::vector<double> r(mNumParameters);

	u[mNumParameters]=0.;
	u[mNumParameters+1]=0.;
	u[mNumParameters+2]=0.;

	for(size_t i=0;i<mNumParameters;++i)
		norm+=u[i]*u[i];
	norm=std::sqrt(norm);
	for(size_t i=0;i<mNumParameters;++i)
		u[i]*=1./norm;

	SetMaxGradientCalls((int) mNumParameters);

	while(true)
	{
		numGradientCalls++;
		 if (numGradientCalls>mMaxGradientCalls)
		 {
             returnValue = eOptimizationReturnAttributes::MAXGRADIENTCALLS;
			 break;
			 //return MAXGRADIENTCALLS;
		 }


		// r set to zero in MatVec product
//		std::cout<<"[MisesWielandt] ATTENTION: has to be D^-1.K not K matrix. \n";
//		r=Au;
		mpCallbackHandlerGrid->Gradient(u,r);
		// input p=r , output r_new=D^-1K r,
		// r of matrix D^-1K
//		y'=D^(-1)*y
		if(mObjectiveType==MAX_EIGENVALUE_OF_PRECOND_MATRIX
				|| mObjectiveType==SPECTRAL_RADIUS_OF_PRECOND_MATRIX
				|| mObjectiveType==CONDITION_NUMBER_OF_PRECOND_MATRIX)
		{
			mpCallbackHandlerGrid->Hessian(r);
		}
		// spectral shift -> does not work
		// y=(D^(-1)K-lambda I)x

		//for spectral radius of M=I-PA
//		for(size_t i=0;i<mNumParameters;++i)
//		{
//			r[i]*=-1;
//			r[i]+=u[i];
//		}
        //prevNorm=norm;
		norm=0.;
		prevLambda=lambda;
		lambda=0.;
		for(size_t i=0;i<mNumParameters;++i)
		{
			lambda+=r[i]*u[i];
			norm+=(r[i]*r[i]); //norm
		}
		norm=std::sqrt(norm);
		for(size_t i=0;i<mNumParameters;++i)
			u[i]=r[i]/norm;

//		std::cout <<"[MisesWielandt] lambda - norm "<<lambda<<" "<<norm<<" " <<"\n";

		if(numGradientCalls!=1)
		{
			if((lambda-prevLambda)>=0&&(lambda-prevLambda)<tol)
			{
                returnValue = eOptimizationReturnAttributes::DELTAOBJECTIVEBETWEENCYCLES;
				break;
			}
			else if((lambda-prevLambda)<0&&(lambda-prevLambda)>-tol)
			{
                returnValue = eOptimizationReturnAttributes::DELTAOBJECTIVEBETWEENCYCLES;
				break;
			}
		}
	}
	// so far condition number equal max eigenvalue
	mObjective=lambda;


	mIsBuild = true;

	if (mVerboseLevel>0)
	{
		std::cout<< " "  << std::endl;
		std::cout<< "[MisesWielandt] Number of gradient calls......... " << numGradientCalls << std::endl;
		std::cout <<"[MisesWielandt] Max eigenvalue .................. "<<lambda << "\n";
		std::cout<< "[MisesWielandt] Active convergence criterion..... " ;
		switch (returnValue)
		{
            case eOptimizationReturnAttributes::MAXFUNCTIONCALLS:
				std::cout<< "Maximum number of function calls reached." << std::endl;
				break;
            case eOptimizationReturnAttributes::MAXGRADIENTCALLS:
				std::cout<< "Maximum number of gradient calls reached." << std::endl;
				break;
            case eOptimizationReturnAttributes::MAXHESSIANCALLS:
				std::cout<< "Maximum number of hessian calls reached." << std::endl;
				break;
            case eOptimizationReturnAttributes::MAXITERATIONS:
				std::cout<< "Maximum number of iterations reached." << std::endl;
				break;
            case eOptimizationReturnAttributes::NORMGRADIENT:
				std::cout<< "lambda_min of preconditioned gradient smaller than prescribed value." << std::endl;
				break;
            case eOptimizationReturnAttributes::DELTAOBJECTIVEBETWEENCYCLES:
				std::cout<< "Decrease in mObjective function between two consecutive cycles is smaller than prescribed value."<<std::endl;
				break;
			default:
				std::cout<< "Unknown convergence criterion." << std::endl;
				break;
		}
		std::cout << std::endl;
	}
    return static_cast<int>(returnValue);
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::MisesWielandt::Info () const
{
    NuTo::Optimizer::InfoBase();
	std::cout<< "AccuracyGradient" << mAccuracyGradient << std::endl;
	std::cout<< "MinDeltaObjBetweenRestarts" << mMinDeltaObjBetweenRestarts << std::endl;
	std::cout<< "MaxGradientCalls" << mMaxGradientCalls << std::endl;
	std::cout<< "MaxHessianCalls" << mMaxHessianCalls << std::endl;
	std::cout<< "MaxIterations" << mMaxIterations << std::endl;
	std::cout<< "ShowSteps" << mShowSteps << std::endl;
}
