#include <iomanip>
#include <iostream>
#include "optimize/Jacobi.h"
#include "base/Timer.h"

//! @brief ... Optimize routine - optimize displacement or error according to input
//! @brief ... Optimize residual if $|r|\inequal 0$
//! @brief ... Optimize displacements if |r|=0, forces have to be added (at place of residuals)
int NuTo::Jacobi::Optimize()
{
    std::vector<double>& v = GetParametersVec();
    std::vector<double>& f = mpCallbackHandlerGrid->GetRightHandSide();

    int returnValue = Optimize(v, f);

    return returnValue;
}

//! @brief ... Optimize routine - optimize solution vector (displacements or error)
//! @brief ... equation system: Kv=f or Ke=r
//! @brief ... $v \leftarrow v - \omega D^{-1} K v + \omega D^{-1} f
//! @param ... v - solution vector, f - right hand sight vector
//! @return ... optimization_return_attribute
int NuTo::Jacobi::Optimize(std::vector<double>& v, std::vector<double>& f)
{
    NuTo::Timer timer(__PRETTY_FUNCTION__, mShowTime);

    double rErrorNorm = 0.;
    int numGradientCalls(0), // number of gradient calls
            curIteration(0); // number of iterations

    eOptimizationReturnAttributes returnValue;


    std::vector<double> vNext(mNumParameters + 3, 0.);
    std::vector<double> rhs = f;

    int precision = 6;
    int width = 10;

    bool converged(false);

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
			throw Exception ( "[Jacobi::Save]File type not implemented." );
		}
	}
	catch ( boost::archive::archive_exception &e )
	{
		std::string s ( std::string ( "[Jacobi::Save]File save exception in boost - " ) +std::string ( e.what() ) );
		std::cout << s << "\n";
		throw Exception ( s );
	}
	catch ( Exception &e )
	{
        throw;
	}
	catch ( std::exception &e )
	{
		throw Exception ( e.what() );
	}
	catch ( ... )
	{
		throw Exception ( "[Matrix::Save]Unhandled exception." );
	}
}

    mpCallbackHandlerGrid->Hessian(rhs);
    // weighting factor (normally at Hessian())
    for (size_t i = 0; i < mNumParameters; ++i)
        rhs[i] *= mOmega;

    while (!converged)
    {
        numGradientCalls++;
        if (numGradientCalls > mMaxGradientCalls)
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw Exception ( "[NuTo::Jacobi::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             oba & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
        curIteration++;
        if (curIteration > mMaxIterations)
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw Exception ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

            if ( tmpString!=GetTypeId() )
                throw Exception ( "[NuTo::Jacobi::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             oxa & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
        // reset gNext to zero in gradient
        mpCallbackHandlerGrid->Gradient(v, vNext);
        // multiply with point diagonal preconditoner
        mpCallbackHandlerGrid->Hessian(vNext);
        rErrorNorm = 0.;
        for (size_t i = 0; i < mNumParameters; ++i)
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw Exception ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

            if ( tmpString!=GetTypeId() )
                throw Exception ( "[NuTo::Jacobi::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

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
            throw Exception ( "[Matrix::Restore]File type not implemented" );
        }
    }
    catch ( Exception &e )
    {
        throw;
    }
    catch ( std::exception &e )
    {
        throw Exception ( e.what() );
    }
    catch ( ... )
    {
        throw Exception ( "[Matrix::Restore]Unhandled exception." );
    }
    return static_cast<int>(returnValue);
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::Jacobi::Info() const
{
    //    NuTo::Optimizer::InfoBase();
    std::cout << "AccuracyGradient" << mAccuracyGradient << std::endl;
    std::cout << "MinDeltaObjBetweenRestarts" << mMinDeltaObjBetweenRestarts << std::endl;
    std::cout << "MaxGradientCalls" << mMaxGradientCalls << std::endl;
    std::cout << "MaxHessianCalls" << mMaxHessianCalls << std::endl;
    std::cout << "MaxIterations" << mMaxIterations << std::endl;
    std::cout << "ShowSteps" << mShowSteps << std::endl;
}
