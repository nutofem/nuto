// $Id $


//! @author Andrea Keßler, ISM
//! @brief ... conjugate gradient method without global matrix matrix-vector product

#include "nuto/optimize/ConjugateGradientGrid.h"
#include "nuto/optimize/MisesWielandt.h"
//#include "nuto/optimize/MultiGrid.h"
#define machine_precision 1e-15
//sqrt machine_precision

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
#define tol 1e-8

int NuTo::ConjugateGradientGrid::Optimize()
{
#ifdef SHOW_TIME
    std::clock_t startOpt,endOpt;
    startOpt=clock();
#endif
	double alpha=0,
		   beta=0,
		   alphaNumerator=0,
		   alphaDenominator=0,
		   norm=0.,
		   betaNumerator=0;

    int numFunctionCalls(0),   // number of function calls
		 numGradientCalls(0),   // number of gradient calls
		 numHessianCalls(0),    // number of hessian calls
		 curIteration(0),       //number of iterations
//		 repeatInitial(10),		// number after which the inital error is recomputed
		 curCycle(0);           //number of iterations without restart


	optimization_return_attributes returnValue;
#pragma acc data  copyin()

	// u  = parameters
	// r  = gradient
	// pr = preconditioned gradient = p r
	// d  = search direction
	// h  = scaled search direction = K d

	std::vector<double> &v=mpCallbackHandlerGrid->GetParameters();
	std::vector<double> &f=mpCallbackHandlerGrid->GetRightHandSide();
	std::vector<double> rInitial(mNumParameters+3);
	std::vector<double> r(mNumParameters+3);
	std::vector<double> pr(mNumParameters+3);
//	std::cout<<"pr[numpar]="<<pr[mNumParameters-1]<<"\n";
	std::vector<double> h(mNumParameters+3);
	std::vector<double> d(mNumParameters+3);

#pragma acc data copy(v),copyin(r,pr,h,d)

	int precision = 6;
	int width = 10;
	bool converged(false);
	double rAccuracyGradientScaled = mAccuracyGradient;
	double rAccuracyGradientSquare =rAccuracyGradientScaled*rAccuracyGradientScaled;
//	if(mVerboseLevel>2)
//		std::cout<<"[ConjugateGradientGrid::Optimize] gradient accuracy "<<rAccuracyGradientScaled <<std::endl;


	int localMaxGradientCalls=(int) mNumParameters/2;
//	int localMaxGradientCalls=2*mNumParameters;
	if (localMaxGradientCalls<mMaxGradientCalls)
		SetMaxGradientCalls(localMaxGradientCalls);

	// calculate objective
	if (numFunctionCalls++>mMaxFunctionCalls)
	{
		converged = true;
		returnValue = MAXFUNCTIONCALLS;
	}


	double residualNorm0=0.;
	while(!converged)
	{
		numGradientCalls++;
		 if (numGradientCalls>mMaxGradientCalls)
		 {
			 converged = true;
			 returnValue = MAXGRADIENTCALLS;
			 break;
		 }


		// all  mNumParameters repeat start equations
//		if (curCycle%repeatInial==0)
		// no upate with start equations
		if (curCycle==0)
		{
			if (numHessianCalls>mMaxHessianCalls)
			{
				converged = true;
				returnValue = MAXHESSIANCALLS;
				break;
			}

			double help=0.0;
			for(size_t i=0;i<mNumParameters;++i)
			{
				help+=v[i]*v[i];
			}
//			std::cout<<__FILE__<<" "<<__LINE__<<" help "<<help<<std::endl;
			// help is only zero when the starting with error equation Ke=r
			if (help!=0)
			{
				// set h to zero done in matrix-vector routine
				mpCallbackHandlerGrid->Gradient(v,r);

				for(size_t i=0;i<mNumParameters;++i)
				{
//				if(help==0) //calculate r=f-Ku, multiply -1, else r=Ke positive Gradient needed
					r[i]*=-1;
				}
			}
			alphaNumerator=0;
			residualNorm0=0;
			for(size_t i=0;i<mNumParameters;++i)
			{
				r[i]+=f[i];
				pr[i]=r[i];
			}
			mpCallbackHandlerGrid->Hessian(pr);
			++numHessianCalls;

			for(size_t i=0;i<mNumParameters;++i)
			{
				//upate for start solution
 				d[i]=pr[i];
				alphaNumerator+=r[i]*pr[i];
				residualNorm0+=r[i]*r[i];
			}

//			residualNorm0=alphaNumerator;
		//	curCycle = 0;
//			std::cout<<" residualNorm0 "<<residualNorm0<<"\n";
		}

		// set h to zero done in matrix-vector routine
		mpCallbackHandlerGrid->Gradient(d,h);

	    alphaDenominator=0.;
//	    squaredNorm=0;
		for(size_t i=0;i<mNumParameters;++i)
			alphaDenominator+=d[i]*h[i];
		alpha = alphaNumerator/alphaDenominator;
		// neg/zero curvature of direction
//		if(alphaDenominator<=0) //www.weizmann
//		{
//			std::cout<< "[ConjugateGradientGrid::Optimize] Negative curacture of search direction. " << std::endl;
//			double normDirection=0.;
//			for(size_t i=0;i<mNumParameters;++i)
//				normDirection+=d[i]*d[i];
//			for(size_t i=0;i<mNumParameters;++i)
//				v[i]=d[i]/normDirection;
//			converged=true;
//			break;
//		}
//
		betaNumerator =0;
		norm =0;
		for(size_t i=0;i<mNumParameters;++i)
		{
			v[i]+= alpha*d[i];
			r[i]-=alpha*h[i];
			pr[i]=r[i];
			norm+=r[i]*r[i];
		}
		mpCallbackHandlerGrid->Hessian(pr);
		++numHessianCalls;

		for(size_t i=0;i<mNumParameters;++i)
			betaNumerator+=r[i]*pr[i];

//		if (mVerboseLevel>0)
//		{
//			std::cout<<"Cycle "<<curCycle<<"\n";
//			std::cout<<"CG r ";
//			for(size_t i=0;i<mNumParameters;++i)
//				std::cout<<r[i]<<" ";
//			std::cout<<"\n";
//
//			std::cout<<"CG pr ";
//			for(size_t i=0;i<mNumParameters;++i)
//				std::cout<<pr[i]<<" ";
//			std::cout<<"\n";
//			std::cout<<"CG v ";
//			for(size_t i=0;i<mNumParameters;++i)
//				std::cout<<v[i]<<" ";
//			std::cout<<"\n";
//		}

		beta = betaNumerator/ alphaNumerator;
		alphaNumerator = betaNumerator;

		if (beta<0)
		{
//			std::cout<< "[ConjugateGradientGrid::Optimize] Set beta ("<< beta <<") to zero not done" << std::endl;
			std::cout<< "[ConjugateGradientGrid::Optimize] Set beta ("<< beta <<") to zero " << std::endl;
			beta=0;
		}

		for(size_t i=0;i<mNumParameters;++i)
		{
			d[i] *=beta;
			d[i] +=pr[i];
		}


		// new convergence criteria
//		if(curCycle>0 && curCycle%repeatInitial==0)
//		{
//			// set h to zero done in matrix-vector routine
//			mpCallbackHandlerGrid->Gradient(v,rInitial);
//
//			for(size_t i=0;i<mNumParameters;++i)
//			{
//				rInitial[i]*=-1;
//				//rInitial[i]+=f[i];
//				pr[i]=rInitial[i];
//			}
////			mpCallbackHandlerGrid->Hessian(pr);
////			std::cout<<"\n";
////			std::cout << " rIy ["<<curCycle<<"] ";
////			for(size_t i=1;i<mNumParameters;i+=3)
////				std::cout<<r[i]<<" ";
//
//			double normInitial=0.;
//			for(size_t i=0;i<mNumParameters;++i)
//			{
//				normInitial+=pr[i]*rInitial[i];
//			}
//			std::cout<<"normInitial "<<normInitial<<"\n";
//			if(2*betaNumerator<normInitial)
//			{
//				converged = true;
//				returnValue = DELTAOBJECTIVEBETWEENCYCLES;
//			}
//
//		}

//		std::cout << "betaNumerator "<<betaNumerator<< " tol: "<<rAccuracyGradientScaled*rAccuracyGradientScaled*residualNorm0<< std::endl;
		// relative norm
//		std::cout<<" relative norm \n";
//		betaNumerator/=residualNorm0;
		norm/=residualNorm0;

//		if ((betaNumerator>0 && betaNumerator<rAccuracyGradientSquare) ||
//				(betaNumerator<0 && betaNumerator>-rAccuracyGradientSquare))
		//increase iteration and curCycle
		curCycle++;
		curIteration++;
		if (norm<rAccuracyGradientSquare)
		{
			converged = true;
			returnValue = NORMGRADIENT;
		}


		if (curIteration>mMaxIterations)
		{
			converged = true;
			returnValue = MAXITERATIONS;
			break;
		}

		numFunctionCalls++;
		if (numFunctionCalls>mMaxFunctionCalls)
		{
			converged = true;
			returnValue = MAXFUNCTIONCALLS;
			break;
		}
	}

	isBuild = true;

#ifdef SHOW_TIME
    endOpt=clock();
    if (mShowTime && mVerboseLevel>0)
    {
    	std::cout<< "[ConjugateGradientGrid::Optimize] Elapsed time (sec)............. " << difftime(endOpt,startOpt)/CLOCKS_PER_SEC << std::endl;
		std::ofstream file;
		file.open("sumOutput",std::ofstream::out|std::ofstream::app);
		if(file)
		{
			//output : voxels in one direction - dofs -
			// nbr grids - solMeth - nbr cycles -nbr pre -nbr post - time -its
			file<< difftime(endOpt,startOpt)/CLOCKS_PER_SEC <<" "<<curIteration<<"\n";
			file.close();
		}
    }
#endif
	if (mVerboseLevel>0)
	{
		std::cout<< " "  << std::endl;
		std::cout<< "[ConjugateGradientGrid::Optimize] "  << std::endl;
		std::cout<< "Number of Iterations............. " << curIteration << std::endl;
		std::cout<< "Active convergence criterion..... " ;
		switch (returnValue)
		{
			case MAXFUNCTIONCALLS:
				std::cout<< "Maximum number of function calls reached." << std::endl;
				break;
			case MAXGRADIENTCALLS:
				std::cout<< "Maximum number of gradient calls reached." << std::endl;
				break;
			case MAXHESSIANCALLS:
				std::cout<< "Maximum number of hessian calls reached." << std::endl;
				break;
			case MAXITERATIONS:
				std::cout<< "Maximum number of iterations reached." << std::endl;
				break;
			case NORMGRADIENT:
				std::cout<< "Norm of preconditioned gradient smaller than prescribed value." << std::endl;
				break;
			case DELTAOBJECTIVEBETWEENCYCLES:
				std::cout<< "Residual discrepancy is smaller." << std::endl;
				break;
			default:
				std::cout<< "Unknown convergence criterion." << std::endl;
				break;
		}
		std::cout << std::endl;

//		std::cout.precision(precision);
//		std::cout << std::setw(width)<< "[ConjugateGradientGrid::Optimize] displacements " ;
//		for (size_t count=0; count<mNumParameters; count++)
//		{
//			std::cout << std::setw(width)<<v[count] << "   " ;
//
//		}
//		std::cout << std::endl;
	}
	objective=sqrt(betaNumerator);
	return returnValue;
}

#ifdef ENABLE_SERIALIZATION
//! @brief ... save the object to a file
//! @param filename ... filename
//! @param rType ... type of file, either BINARY, XML or TEXT
void NuTo::ConjugateGradientGrid::Save ( const std::string &filename, std::string rType)const
{
	try
	{
		//transform to uppercase
		std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
		std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
		std::string tmpStr ( GetTypeId() );
		std::string baseClassStr = tmpStr.substr ( 4,100 );
		if (rType=="BINARY")
		{
			boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
			oba & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oba & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                & BOOST_SERIALIZATION_NVP(mMaxIterations)
                & BOOST_SERIALIZATION_NVP(mShowSteps);
		}
		else if (rType=="XML")
		{
			boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
			oxa & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oxa & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                & BOOST_SERIALIZATION_NVP(mMaxIterations)
                & BOOST_SERIALIZATION_NVP(mShowSteps);
		}
		else if (rType=="TEXT")
		{
			boost::archive::text_oarchive ota ( ofs, std::ios::binary );
			ota & boost::serialization::make_nvp ( "Object_type", tmpStr );
            ota & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                & BOOST_SERIALIZATION_NVP(mMaxIterations)
                & BOOST_SERIALIZATION_NVP(mShowSteps);
		}
		else
		{
			throw MathException ( "[ConjugateGradientGrid::Save]File type not implemented." );
		}
	}
	catch ( boost::archive::archive_exception &e )
	{
		std::string s ( std::string ( "[ConjugateGradientGrid::Save]File save exception in boost - " ) +std::string ( e.what() ) );
		std::cout << s << "\n";
		throw MathException ( s );
	}
	catch ( MathException &e )
	{
		throw e;
	}
	catch ( std::exception &e )
	{
		throw MathException ( e.what() );
	}
	catch ( ... )
	{
		throw MathException ( "[Matrix::Save]Unhandled exception." );
	}
}


//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::ConjugateGradientGrid::Restore ( const std::string &filename,  std::string rType)
{
    try
    {
		//transform to uppercase
		std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        std::string tmpString;
		if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw OptimizeException ( "[NuTo::ConjugateGradientGrid::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             oba & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
		else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MathException ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

            if ( tmpString!=GetTypeId() )
                throw OptimizeException ( "[NuTo::ConjugateGradientGrid::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             oxa & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
		else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MathException ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

            if ( tmpString!=GetTypeId() )
                throw OptimizeException ( "[NuTo::ConjugateGradientGrid::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             ota & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
		else
		{
            throw MathException ( "[Matrix::Restore]File type not implemented" );
        }
    }
    catch ( MathException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MathException ( e.what() );
    }
    catch ( ... )
    {
        throw MathException ( "[Matrix::Restore]Unhandled exception." );
    }
}
#endif // ENABLE_SERIALIZATION

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name ConjugateGradientGrid
std::string NuTo::ConjugateGradientGrid::GetTypeId()const
{
    return std::string("ConjugateGradientGrid");
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::ConjugateGradientGrid::Info () const
{
	std::cout<< "ConjugateGradientGrid\n";
	std::cout<<"-----------------------------------------------------------------------------------\n";
	std::cout<< "AccuracyGradient ....................... " << mAccuracyGradient << std::endl;
	std::cout<< "MaxGradientCalls ....................... " << mMaxGradientCalls << std::endl;
	std::cout<< "MaxHessianCalls ........................ " << mMaxHessianCalls << std::endl;
	std::cout<< "MaxIterations .......................... " << mMaxIterations << std::endl;
}
