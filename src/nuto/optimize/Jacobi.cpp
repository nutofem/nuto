// $Id $


//! @author Andrea Keßler, ISM
//! @brief ...jacobi method


#include "nuto/optimize/Jacobi.h"
#include "nuto/math/FullMatrix.h"
#define machine_precision 1e-15
//sqrt machine_precision

#define tol 1e-8

//! @brief ... Optimize routine - optimize displacement or residual according to input
//! @brief ... Optimize residual if $|r|\inequal 0$
//! @brief ... Optimize displacements if |r|=0
int NuTo::Jacobi::Optimize()
{

#ifdef SHOW_TIME
    std::clock_t startJac,endJac;
    startJac=clock();
#endif
//		std::cout<< __FILE__<<" "<<__LINE__<< "[Jacobi::Optimize]" << std::endl;
    double rErrorNorm=0.;
    int  numGradientCalls(0),   // number of gradient calls
		 curIteration(0);      //number of iterations


	optimization_return_attributes returnValue;

	std::fstream outputTime;
	std::string filename = "timeOutput";

	// u  = parameters
	// p  = preconditioner D^-1 as vector

	std::vector<double> &u=mpCallbackHandlerGrid->GetParameters();
	std::vector<double> &r=mpCallbackHandlerGrid->GetResidual();
	std::vector<double> gNext(mNumParameters,0.);
	std::vector<double> g(mNumParameters,0.);
	std::vector<double> p(mNumParameters,0.);

	int precision = 6;
	int width = 10;

	bool converged(false);

	int localMaxGradientCalls=(int) mNumParameters;
//	int localMaxGradientCalls=2*mNumParameters;
	if (localMaxGradientCalls<mMaxGradientCalls)
		SetMaxGradientCalls(localMaxGradientCalls);



	boost::dynamic_bitset<> rDofIsConst=mpCallbackHandlerGrid->GetDisplacementConstaints();

	double help=0;
	for(size_t i=0;i<mNumParameters;++i)
	{
		help+=r[i]*r[i];
	}
	if (help!=0) //residual vector is not zero, optimize residual
	{
		for(size_t i=0;i<mNumParameters;++i)
			g[i]=r[i];
	}
	else //optimize displacements
	{
		for(size_t i=0;i<mNumParameters;++i)
			g[i]=u[i];
	}



	while(!converged)
	{
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
		 // reset gNext to zero in gradient
		 mpCallbackHandlerGrid->Gradient(g,gNext);
		// multiply with point diagonal preconditoner
		mpCallbackHandlerGrid->Hessian(gNext);
		rErrorNorm=0.;
		for(size_t i=0;i<mNumParameters;++i)
		{
//			if(!rDofIsConst[i])
			// scaling direct on preconditioner
			// plus omega?
			// do not forget sign!!!
			gNext[i]*=-mOmega;
			gNext[i]+=g[i];
			// if(help==0)
			// gNext[i]+=mOmega*p[i]*f[i]; for forces
			rErrorNorm+=(gNext[i]-g[i])*(gNext[i]-g[i]);
			g[i]=gNext[i];
		}

		if (rErrorNorm < mAccuracyGradient*mAccuracyGradient)
		{
			converged = true;
			returnValue = DELTAOBJECTIVEBETWEENCYCLES;
			break;
		}

//		if (mVerboseLevel>4 && curIteration==100)
//		{
//			std::cout.precision(precision);
//			std::cout <<std::setw(width)<<"[Jacobi::Optimize]  It.: "<< curIteration<< " norm squared "<<rErrorNorm<<std::endl;
//		if (help!=0)
//			std::cout << std::setw(width)<< "[Jacobi::Optimize] residual " ;
//		else
//			std::cout << std::setw(width)<< "[Jacobi::Optimize] displacements " ;
//			for (size_t count=0; count<mNumParameters; count++)
//			{
//				std::cout << std::setw(width)<< g[count] << "   " ;
//			}
//			std::cout << std::endl;
//		}
	}

	isBuild = true;
	if (help!=0) //residual vector is not zero, optimize residual
	{
		for(size_t i=0;i<mNumParameters;++i)
			r[i]=g[i];
	}
	else //optimize displacements
	{
		for(size_t i=0;i<mNumParameters;++i)
			u[i]=g[i];
	}

#ifdef SHOW_TIME
    endJac=clock();
	if (mShowTime &&mVerboseLevel>0 )
		std::cout<< "[Jacobi::Optimize] Elapsed time (sec)............. " << difftime(endJac,startJac)/CLOCKS_PER_SEC << std::endl;
    outputTime.open(filename,std::fstream::out|std::fstream::app);
 	outputTime<<(difftime(endJac,startJac)/CLOCKS_PER_SEC)<<"   "<<curIteration<<"\n";
	outputTime.close();

#endif
	if (mVerboseLevel>0)
	{
		std::cout<< "[Jacobi::Optimize] Number of Iterations............. " << curIteration << std::endl;
		std::cout<< "[Jacobi::Optimize] Active convergence criterion..... " ;
		switch (returnValue)
		{
			case MAXGRADIENTCALLS:
				std::cout<< "Maximum number of gradient calls reached." << std::endl;
				break;
			case MAXITERATIONS:
				std::cout<< "Maximum number of iterations reached." << std::endl;
				break;
			case DELTAOBJECTIVEBETWEENCYCLES:
				std::cout<< "Norm of error smaller than prescribed value." << std::endl;
				objective=sqrt(rErrorNorm);
				break;
			default:
				std::cout<< "Unknown convergence criterion." << std::endl;
				break;
		}
		std::cout << std::endl;

		std::cout.precision(precision);
		if(help==0)
			std::cout << std::setw(width)<< "[Jacobi::Optimize] displacements optimized\n "
					"displacements " ;
		else
			std::cout << std::setw(width)<< "[Jacobi::Optimize] residual optimized\n "
					"residual " ;
		for (size_t count=0; count<mNumParameters; count++)
		{
			std::cout << std::setw(width)<<g[count] << "   " ;

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
//    NuTo::Optimizer::InfoBase();
	std::cout<< "AccuracyGradient" << mAccuracyGradient << std::endl;
	std::cout<< "MinDeltaObjBetweenRestarts" << mMinDeltaObjBetweenRestarts << std::endl;
	std::cout<< "MaxGradientCalls" << mMaxGradientCalls << std::endl;
	std::cout<< "MaxHessianCalls" << mMaxHessianCalls << std::endl;
	std::cout<< "MaxIterations" << mMaxIterations << std::endl;
	std::cout<< "ShowSteps" << mShowSteps << std::endl;

}
