// $Id: RungeKuttaDormandPrince.cpp 575 2011-09-20 18:05:35Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#endif // ENABLE_SERIALIZATION

# ifdef _OPENMP
#include <omp.h>
# endif

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/RungeKuttaDormandPrince.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKuttaDormandPrince::RungeKuttaDormandPrince (StructureBase* rStructure)  : RungeKuttaBase (rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKuttaDormandPrince::Info()const
{
	TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines (this is wrong do)
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKuttaDormandPrince::CalculateCriticalTimeStep()const
{
	double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
	return 2.8/std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKuttaDormandPrince::GetStageTimeFactor(int rStage)const
{
	assert(rStage<7);
	double s;
	switch(rStage)
	{
	case 0:
		s = 0.;
		break;
	case 1:
		s = 0.2;
		break;
	case 2:
		s = 0.3;
		break;
	case 3:
		s = 0.8;
		break;
	case 4:
		s = 8./9.;
		break;
	case 5:
		s = 1.0;
		break;
	case 6:
		s = 1.0;
		break;
	default:
		throw MechanicsException("[NuTo::RungeKuttaDormandPrince::GetStageTimeFactor] rStage<7.");
	}
	return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKuttaDormandPrince::HasTimeChanged(int rStage)const
{
	assert(rStage<7);
	bool s;
	switch(rStage)
	{
	case 0:
		s = false; //same as last step from the last iteration
		break;
	case 1:
		s = true;
		break;
	case 2:
		s = true;
		break;
	case 3:
		s = true;
		break;
	case 4:
		s = true;
		break;
	case 5:
		s = true;
		break;
	case 6:
		s = false;
		break;
	default:
		throw MechanicsException("[NuTo::RungeKuttaDormandPrince::HasTimeChanged] rStage<7.");
	}
	return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKuttaDormandPrince::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage)const
{
	assert(rStage<7);
	assert(rWeight.size()==6);
	switch(rStage)
	{
	case 0:
		break;
	case 1:
		rWeight[0] = 0.2;
		break;
	case 2:
		rWeight[0] = 3./40.;
		rWeight[1] = 9./40.;
		break;
	case 3:
		rWeight[0] = 44./45.;
		rWeight[1] = -56./15.;
		rWeight[2] = 32./9.;
		break;
	case 4:
		rWeight[0] = 19372./6561.;
		rWeight[1] = -25360./2187.;
		rWeight[2] = 64448./6561.;
		rWeight[3] = -212./729.;
		break;
	case 5:
		rWeight[0] = 9017./3168.;
		rWeight[1] = -355./33.;
		rWeight[2] = 46732./5247.;
		rWeight[3] = 49./176.;
		rWeight[4] = -5103./18656.;
		break;
	case 6:
		rWeight[0] = 35./384.;
		rWeight[1] = 0.;
		rWeight[2] = 500./1113.;
		rWeight[3] = 125./192.;
		rWeight[4] = -2187./6784.;
		rWeight[5] = 11./84.;
		break;
	default:
		throw MechanicsException("[NuTo::RungeKuttaDormandPrince::GetStageDerivativeFactor] rStage<7.");
	}
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKuttaDormandPrince::GetStageWeights(int rStage)const
{
	assert(rStage<7);
	double s;
	switch(rStage)
	{
	case 0:
		s = 35./385;
		break;
	case 1:
		s = 0.;
		break;
	case 2:
		s = 500./1113.;
		break;
	case 3:
		s = 125./192.;
		break;
	case 4:
		s = -2187./6784.;
		break;
	case 5:
		s = 11./84.;
		break;
	case 6:
		s = 0.;
		break;
	default:
		throw MechanicsException("[NuTo::RungeKuttaDormandPrince::GetStageWeights] rStage<7.");
	}
	return s;
}
#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::RungeKuttaDormandPrince::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaDormandPrince::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaDormandPrince::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaDormandPrince::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaDormandPrince::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaDormandPrince::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::RungeKuttaDormandPrince::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of RungeKuttaDormandPrince" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RungeKuttaBase);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of RungeKuttaDormandPrince" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::RungeKuttaDormandPrince::GetTypeId()const
{
    return std::string("RungeKuttaDormandPrince");
}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::RungeKuttaDormandPrince::Restore (const std::string &filename, std::string rType )
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
                throw MechanicsException ( "[RungeKuttaDormandPrince::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[RungeKuttaDormandPrince::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[RungeKuttaDormandPrince::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            ota & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else
        {
            throw MathException ( "[Matrix::Restore]File type not implemented" );
        }
    }
    catch ( MechanicsException &e )
    {
        throw;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
    catch ( ... )
    {
        throw MechanicsException ( "[RungeKuttaDormandPrince::Restore]Unhandled exception." );
    }
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::RungeKuttaDormandPrince::Save (const std::string &filename, std::string rType )const
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
            oba & boost::serialization::make_nvp(tmpStr.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oxa & boost::serialization::make_nvp(tmpStr.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpStr );
            ota & boost::serialization::make_nvp(tmpStr.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[RungeKuttaDormandPrince::Save]File type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[RungeKuttaDormandPrince::Save]File save exception in boost - " ) +std::string ( e.what() ) );
        std::cout << s << "\n";
        throw MathException ( s );
    }
    catch ( MechanicsException &e )
    {
        throw;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
    catch ( ... )
    {
        throw MechanicsException ( "[RungeKuttaDormandPrince::Save] Unhandled exception." );
    }
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::RungeKuttaDormandPrince)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
