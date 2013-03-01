// $Id $


//! @author Andrea Keßler, ISM
//! @brief ...jacobi method


#include "nuto/optimize/Jacobi.h"
#include "nuto/math/FullMatrix.h"
#define machine_precision 1e-15
//sqrt machine_precision

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
#define tol 1e-8

int NuTo::Jacobi::Optimize()
{
#ifdef SHOW_TIME
    std::clock_t startJac,endJac;
    startJac=clock();
#endif
    double rErrorNorm=0.;
    double rResidualNorm=0.;
    int numFunctionCalls(0),   // number of function calls
		 numGradientCalls(0),   // number of gradient calls
		 numHessianCalls(0),    // number of hessian calls
		 curIteration(0),       //number of iterations
		 curCycle(0);           //number of iterations without restart


	optimization_return_attributes returnValue;

	std::fstream outputTime;
	std::string filename = "timeOutput";

	// u  = parameters
	// r  = parameters of next step
	// rprev  = parameters of previous step
	// p  = preconditioner D^-1 as vector

	std::vector<double> &u=mpCallbackHandlerGrid->GetParameters();
	std::vector<double> &r=mpCallbackHandlerGrid->GetResidual();
//	GetParameters(u);
	std::vector<double> uNext(mNumParameters,0.);
	std::vector<double> p(mNumParameters,0.);

	int precision = 6;
	int width = 10;

	if (mVerboseLevel>2)
		std::cout<< __FILE__<<" "<<__LINE__<< "[NuTo::Jacobi::Optimize] numParameters "<< mNumParameters << std::endl;
	bool converged(false);
	std::cout<<"[NuTo::Jacobi::Optimize] gradient accuracy "<<mAccuracyGradient <<std::endl;

	int localMaxGradientCalls=(int) mNumParameters;
//	int localMaxGradientCalls=2*mNumParameters;
	if (localMaxGradientCalls<mMaxGradientCalls)
		SetMaxGradientCalls(localMaxGradientCalls);


	// calculate Diag preonditioner;
	mpCallbackHandlerGrid->Hessian(p);

	boost::dynamic_bitset<> rDofIsConst=mpCallbackHandlerGrid->GetDisplacementConstaints();



	while(!converged)
	{
		if (numFunctionCalls++ > mMaxFunctionCalls)
		{
			converged = true;
			returnValue = MAXFUNCTIONCALLS;
			break;
		}

		numGradientCalls++;
		 if (numGradientCalls>mMaxGradientCalls)
		 {
			 converged = true;
			 returnValue = MAXGRADIENTCALLS;
			 break;
			 //return MAXGRADIENTCALLS;
		 }
		 curIteration++;
		 if (curIteration>mMaxIterations)
		 {
			 converged = true;
			 returnValue = MAXITERATIONS;
			 break;
			 //return MAXGRADIENTCALLS;
		 }


		//calculate gradient
		//force has to be added
		// print r
//		std::cout<<" u ";
//		for(size_t i=0;i<mNumParameters;++i)
//			std::cout<<u[i]<<" ";
//		std::cout<<"\n";



		mpCallbackHandlerGrid->Gradient(u,uNext);
		rResidualNorm=0.;
		rErrorNorm=0.;
		for(size_t i=0;i<mNumParameters;++i)
		{
//			if(!rDofIsConst[i])
////				rResidualNorm+=(-uNext[i])*(-uNext[i]);
//				rResidualNorm+=(r[i]-uNext[i])*(r[i]-uNext[i]);
			uNext[i]*=-mOmega*p[i];
//			if(!rDofIsConst[i])
			uNext[i]+=!rDofIsConst[i]*u[i];
			// if r exist
			uNext[i]+=mOmega*p[i]*r[i];
			if(!rDofIsConst[i])
			{
				rErrorNorm+=(uNext[i]-u[i])*(uNext[i]-u[i]);
				u[i]=uNext[i];
			}
		}
//		std::cout<<"[Jacobi] rErrorNorm "<<rErrorNorm<<" rResidualNorm "<<rResidualNorm<<"\n";
		if (rErrorNorm < mAccuracyGradient*mAccuracyGradient)
		{
			converged = true;
			returnValue = NORMGRADIENT;
			break;
		}

//		if (rResidualNorm < mAccuracyGradient*mAccuracyGradient)
//		{
//			converged = true;
//			returnValue = NORMGRADIENT;
//			break;
//		}

		if (mVerboseLevel>4 && curIteration==100)
		{
			std::cout.precision(precision);
			std::cout <<std::setw(width)<<" It.: "<< curIteration<< " norm squared "<<rErrorNorm<<std::endl;
			std::cout << std::setw(width)<< "u " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< u[count] << "   " ;
			}
			std::cout << std::endl;
		}
	}

	mpCallbackHandlerGrid->Gradient(u,r);
	for(size_t i=0;i<mNumParameters;++i)
		r[i]*=-1;


//	SetParameters(u);
//	mpCallbackHandlerGrid->SetResidual(u);

	isBuild = true;
	objective=sqrt(rErrorNorm);
	std::cout <<" objective " <<objective  <<"\n";

	numFunctionCalls++;

#ifdef SHOW_TIME
    endJac=clock();
    if (mShowTime)
		std::cout<< "Elapsed time (sec)............. " << difftime(endJac,startJac)/CLOCKS_PER_SEC << std::endl;
    outputTime.open(filename,std::fstream::out|std::fstream::app);
 	outputTime<<(difftime(endJac,startJac)/CLOCKS_PER_SEC)<<"   "<<curIteration<<"\n";
	outputTime.close();

#endif
//	if (mVerboseLevel)
//	{
		std::cout<< " "  << std::endl;
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
			default:
				std::cout<< "Unknown convergence criterion." << std::endl;
				break;
		}
		std::cout << std::endl;
//	}
	if(mVerboseLevel>4)
	{
		std::cout.precision(precision);
		std::cout << std::setw(width)<< "displacements " ;
		for (size_t count=0; count<mNumParameters; count++)
		{
			std::cout << std::setw(width)<<u[count] << "   " ;

		}
		std::cout << std::endl;
	}
	return returnValue;
}

#ifdef ENABLE_SERIALIZATION
//! @brief ... save the object to a file
//! @param filename ... filename
//! @param rType ... type of file, either BINARY, XML or TEXT
void NuTo::Jacobi::Save ( const std::string &filename, std::string rType)const
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
			throw MathException ( "[Jacobi::Save]File type not implemented." );
		}
	}
	catch ( boost::archive::archive_exception &e )
	{
		std::string s ( std::string ( "[Jacobi::Save]File save exception in boost - " ) +std::string ( e.what() ) );
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
void NuTo::Jacobi::Restore ( const std::string &filename,  std::string rType)
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
                throw OptimizeException ( "[NuTo::Jacobi::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

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
                throw OptimizeException ( "[NuTo::Jacobi::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

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
                throw OptimizeException ( "[NuTo::Jacobi::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

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
//! @return    class name Jacobi
std::string NuTo::Jacobi::GetTypeId()const
{
    return std::string("Jacobi");
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::Jacobi::Info () const
{
    NuTo::Optimizer::InfoBase();
	std::cout<< "AccuracyGradient" << mAccuracyGradient << std::endl;
	std::cout<< "MinDeltaObjBetweenRestarts" << mMinDeltaObjBetweenRestarts << std::endl;
	std::cout<< "MaxGradientCalls" << mMaxGradientCalls << std::endl;
	std::cout<< "MaxHessianCalls" << mMaxHessianCalls << std::endl;
	std::cout<< "MaxIterations" << mMaxIterations << std::endl;
	std::cout<< "ShowSteps" << mShowSteps << std::endl;
}
