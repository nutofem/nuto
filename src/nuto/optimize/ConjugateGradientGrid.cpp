// $Id $


//! @author Andrea Keßler, ISM
//! @brief ... conjugate gradient method without global matrix matrix-vector product


#include "nuto/optimize/ConjugateGradientGrid.h"
#include "nuto/optimize/MisesWielandt.h"
#include "nuto/optimize/MultiGrid.h"
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
		   betaNumerator=0;
//		   squaredNorm=0;

    int numFunctionCalls(0),   // number of function calls
		 numGradientCalls(0),   // number of gradient calls
		 numHessianCalls(0),    // number of hessian calls
		 curIteration(0),       //number of iterations
		 curCycle(0);           //number of iterations without restart


	optimization_return_attributes returnValue;

	std::fstream outputTime;
	std::string filename = "timeOutput";

	// u  = parameters
	// r  = gradient
	// pr = preconditioned gradient = p r
	// d  = search direction
	// h  = scaled search direction = K d

	std::vector<double> &u=mpCallbackHandlerGrid->GetParameters();
	std::vector<double> &r=mpCallbackHandlerGrid->GetResidual();
	std::vector<double> pr(mNumParameters);
	std::vector<double> p(mNumParameters,1);
	std::vector<double> h(mNumParameters);
	std::vector<double> d(mNumParameters+3);


	int precision = 6;
	int width = 10;

	if (mVerboseLevel>2)
		std::cout<< __FILE__<<" "<<__LINE__<< "[NuTo::ConjugateGradientGrid::Optimize] numParameters "<< mNumParameters << std::endl;
	bool converged(false);
	double rAccuracyGradientScaled = mAccuracyGradient;
	std::cout<<"[NuTo::ConjugateGradientGrid::Optimize] gradient accuracy "<<rAccuracyGradientScaled <<std::endl;

	int localMaxGradientCalls=(int) mNumParameters;
//	int localMaxGradientCalls=2*mNumParameters;
	if (localMaxGradientCalls<mMaxGradientCalls)
		SetMaxGradientCalls(localMaxGradientCalls);

	// calculate objective
	numFunctionCalls++;
	if (numFunctionCalls>mMaxFunctionCalls)
	{
		converged = true;
		returnValue = MAXFUNCTIONCALLS;
	}


	//0 set if no precondition
//	outputTime.open(filename,std::fstream::out|std::fstream::app);
//	outputTime<<" 0 ";
//	outputTime.close();
	// calculate Diag preonditioner;
	mpCallbackHandlerGrid->Hessian(p);

	//1 set if no scaling
//	std::fstream outputTime;
//	std::string filename = "timeOutput";
//	outputTime.open(filename,std::fstream::out|std::fstream::app);
//	outputTime<<" 1 ";
//	outputTime.close();
	CalcScalingFactors(numHessianCalls,p);

//	 print p
//	std::cout<<" p ";
//	for(int i=0;i<mNumParameters;i+=3)
//		std::cout<<p[i]<<" ";
//	std::cout<<"\n";


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

//determine gradient

		// all  mNumParameters repeat start equations
		//if (curCycle%mNumParameters==0)
		// no upate with start equations
		if (curCycle==0)
		{
			if (numHessianCalls>mMaxHessianCalls)
			{
				converged = true;
				returnValue = MAXHESSIANCALLS;
				break;
			}

			//calculate gradient as a start solution
			if (mVerboseLevel>3)
				std::cout<<"[NuTo::ConjugateGradientGrid::Optimize] ("<<__LINE__<<") calc start direction"<<std::endl;
			double help=0.0;
			for(size_t i=0;i<mNumParameters;++i)
			{
				help+=r[i]*r[i];
			}
//			std::cout<<__FILE__<<" "<<__LINE__<<" help "<<help<<std::endl;
			if (help==0)
				mpCallbackHandlerGrid->Gradient(u,r);
			else
				std::cout<<"[NuTo::ConjugateGradientGrid::Optimize] ("<<__LINE__<<") no gradient calcuation"<<std::endl;

			// have to change sign --> next loop
			// add external load

			// print r
//			std::cout<<" r ";
//			for(size_t i=0;i<mNumParameters;++i)
//				std::cout<<r[i]<<" ";
//			std::cout<<"\n";

//			pr=r;
			for(size_t i=0;i<mNumParameters;++i)
			{
				r[i]*=-1;
//				pr[i]=r[i];
 				pr[i]=p[i]*r[i];
				//upate for start solution
 				d[i]=pr[i];
			}

			alphaNumerator=0;
			for(size_t i=0;i<mNumParameters;++i)
 				alphaNumerator+=r[i]*pr[i];

			if (mVerboseLevel>2)
				std::cout<<__FILE__<<" "<<__LINE__<<" normGradient "<<alphaNumerator << " accuracy " <<rAccuracyGradientScaled<<std::endl;

			if (mVerboseLevel>5 && curCycle>0)
				std::cout<< "   Restart after " <<curCycle << " cycles " << std::endl;
			curCycle = 0;
		std::cout.precision(precision);
//		std::cout <<std::setw(width)<<"It.: "<<curIteration<<" i=3 -> r "<<r[3]<<" p "<<p[3]<<" pr "<<pr[3]<<" a1 "<<alphaNumerator<<"\n";
		}


		//begin with cycles bigger 0

		if (mVerboseLevel>4)
			std::cout<<__FILE__<<" "<<__LINE__<<" calc search direction"<<std::endl;

		// set h to zero done in matrix-vector routine
		mpCallbackHandlerGrid->Gradient(d,h);

	    alphaDenominator=0.;
//	    squaredNorm=0;
		for(size_t i=0;i<mNumParameters;++i)
			alphaDenominator+=d[i]*h[i];
		alpha = alphaNumerator/alphaDenominator;
		betaNumerator =0;
		for(size_t i=0;i<mNumParameters;++i)
		{
			u[i]+= alpha*d[i];
			r[i]-=alpha*h[i];
			pr[i]=p[i]*r[i];
			betaNumerator+=r[i]*pr[i];
		}
		beta = betaNumerator/ alphaNumerator;
		alphaNumerator = betaNumerator;
		for(size_t i=0;i<mNumParameters;++i)
		{
			d[i] *=beta;
			d[i] +=pr[i];
		}

//		std::cout <<std::setw(width)<<"It.: "<<curIteration<<" i=3 -> h "<<h[3]<<" alpha "<<alpha<<" r "<<r[3]<<" p "<<p[3]<<" pr "<<pr[3]<<" b1 "<<betaNumerator<<" d "<<d[3]<<" u "<<u[3]<<"\n";


		if (mVerboseLevel>4)
		{
			std::cout.precision(precision);
			std::cout <<std::setw(width)<<" It.: "<< curIteration<< " norm gradient squared (betaNumerator)"<<betaNumerator<<std::endl;
			std::cout << std::setw(width)<< "pr " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< pr[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "r " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< r[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "h " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< h[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "d " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< d[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "u " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<<u[count] << "   " ;
			}
			std::cout << std::endl;

			std::cout << std::setw(width)<< "alpha "<<alpha<< "beta "<<beta << std::endl;
		}

		if ((betaNumerator>0 && betaNumerator<rAccuracyGradientScaled*rAccuracyGradientScaled) ||
				(betaNumerator<0 && betaNumerator>-rAccuracyGradientScaled*rAccuracyGradientScaled))
		{
			if (mVerboseLevel>2)
			{
				std::cout<< "CONVERGED " << std::endl;
				std::cout<< "Iteration " << curIteration <<" with norm grad squared" << betaNumerator << std::endl;
			}
			converged = true;
			returnValue = NORMGRADIENT;
			break;
		}

		if (beta<0)
		{
			std::cout<< "Set beta ("<< beta <<") to zero " << std::endl;
			beta=0;
		}


//		if (mVerboseLevel>1 && curIteration%mShowSteps==0)
//		if ( curIteration%mShowSteps==0)
//		{
//			std::cout<< "Iteration " << curIteration <<" with norm grad squared" << betaNumerator << std::endl;
//		}

//		std::cout<<" u ";
//		for(size_t i=0;i<mNumParameters;i+=3)
//			std::cout<<u[i]<<" ";
//		std::cout<<"\n";

		//increase iteration and curCycle
		curCycle++;
		curIteration++;

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

	// Get force vector
//	r.clear();
//	CalculateReactionForcesEBE(u,r);
//	std::cout<<" Reaction forces \n";
//	for (size_t count=0; count<mNumParameters; count++)
//	{
//		std::cout << std::setw(width)<< r[count] << "\n   " ;
//	}
//	std::cout<<"\n";

#ifdef SHOW_TIME
    endOpt=clock();
    if (mShowTime)
		std::cout<< "Elapsed time (sec)............. " << difftime(endOpt,startOpt)/CLOCKS_PER_SEC << std::endl;
    outputTime.open(filename,std::fstream::out|std::fstream::app);
 	outputTime<<(difftime(endOpt,startOpt)/CLOCKS_PER_SEC)<<"   "<<curIteration<<"\n";
	outputTime.close();

#endif
//	if (mVerboseLevel)
//	{
		std::cout<< " "  << std::endl;
		std::cout<< "Number of Iterations............. " << curIteration << std::endl;
//		std::cout<< "Number of Function Calls......... " << numFunctionCalls << std::endl;
//		std::cout<< "Number of Gradient Calls......... " << numGradientCalls << std::endl;
//		std::cout<< "Number of Hessian Calls.......... " << numHessianCalls << std::endl;
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
	objective=sqrt(betaNumerator);
	return returnValue;
}

void NuTo::ConjugateGradientGrid::CalcScalingFactors(int& numHessianCalls,std::vector<double> &p)
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif
    //diagonal scaling with scaling factor
	++numHessianCalls;
//    double scalefactor=1;
    double scalefactor=0.0000000001;
	if (mUseMisesWielandt)
	{
		if(!(mpCallbackHandlerGrid->GetWeightingFactor()))
		{
			// scale factor is 2/(2-lambda_max-lambda_min) [Meister: Num. lin. GLS] for Jacobi-Relaxation-Method
			NuTo::MisesWielandt myEigenCalculator(mNumParameters);
			myEigenCalculator.SetVerboseLevel(5);
			myEigenCalculator.SetCallback((mpCallbackHandlerGrid));
			myEigenCalculator.Optimize();
			double lambda_max=myEigenCalculator.GetObjective();
//			double lambda_min=lambda_max/2.;
			// Jacobi-Relaxation-weighting
	//		scalefactor=2./(2-lambda_max-lambda_min);
			// my scale factor
	//		scalefactor=2./(lambda_max*lambda_max);
			scalefactor=1./lambda_max;
			// damping Jacobi: lampda of D-1 K Arbenz_2007
	//		scalefactor=4./(3.*lambda_max);

			mpCallbackHandlerGrid->SetWeightingFactor(scalefactor);
		}
		else
			scalefactor=mpCallbackHandlerGrid->GetWeightingFactor();
	}

	std::fstream outputTime;
	std::string filename = "timeOutput";
    outputTime.open(filename,std::fstream::out|std::fstream::app);
    outputTime<<scalefactor<<"  ";
    outputTime.close();
//    if(mVerboseLevel>0)
	std::cout<<"[ConjugateGradientGrid::CalcScalingFactors] scale factor "<<scalefactor<<"\n";
    for (size_t count=0; count<mNumParameters; ++count)
        p[count] *=scalefactor;

//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::ConjugateGradientGrid::CalcScalingFactors] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
//#endif
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
    NuTo::Optimizer::InfoBase();
	std::cout<< "AccuracyGradient" << mAccuracyGradient << std::endl;
	std::cout<< "MinDeltaObjBetweenRestarts" << mMinDeltaObjBetweenRestarts << std::endl;
	std::cout<< "MaxGradientCalls" << mMaxGradientCalls << std::endl;
	std::cout<< "MaxHessianCalls" << mMaxHessianCalls << std::endl;
	std::cout<< "MaxIterations" << mMaxIterations << std::endl;
	std::cout<< "ShowSteps" << mShowSteps << std::endl;
}
