// $Id: RungeKutta4.cpp 575 2011-09-20 18:05:35Z unger3 $

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
#include "nuto/mechanics/timeIntegration/RungeKutta4.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKutta4::RungeKutta4 (StructureBase* rStructure)  : RungeKuttaBase (rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKutta4::Info()const
{
	TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKutta4::CalculateCriticalTimeStep()const
{
	double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
	return 2.8/std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKutta4::GetStageTimeFactor(int rStage)const
{
	assert(rStage<4);
	double s;
	switch(rStage)
	{
	case 0:
		s = 0.;
		break;
	case 1:
		s = 0.5;
		break;
	case 2:
		s = 0.5;
		break;
	case 3:
		s = 1.0;
		break;
	default:
        throw MechanicsException ( "[NuTo::RungeKutta4::GetStageTimeFactor] rStage>3 not implemented." );
	}
	return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKutta4::HasTimeChanged(int rStage)const
{
	assert(rStage<4);
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
		s = false;
		break;
	case 3:
		s = true;
		break;
	default:
        throw MechanicsException ( "[NuTo::RungeKutta4::HasTimeChanged] rStage>3 not implemented." );
	}
	return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKutta4::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage)const
{
	assert(rStage<4);
	assert(rWeight.size()==3);
	switch(rStage)
	{
	case 0:
		break;
	case 1:
		rWeight[0] = 0.5;
		break;
	case 2:
		rWeight[0] = 0.0;
		rWeight[1] = 0.5;
		break;
	case 3:
		rWeight[0] = 0.0;
		rWeight[1] = 0.0;
		rWeight[2] = 1.0;
		break;
	default:
        throw MechanicsException ( "[NuTo::RungeKutta4::GetStageDerivativeFactor] rStage>3 not implemented." );
	}
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKutta4::GetStageWeights(int rStage)const
{
	assert(rStage<4);
	double s;
	switch(rStage)
	{
	case 0:
		s = 1./6.;
		break;
	case 1:
		s = 1./3.;
		break;
	case 2:
		s = 1./3.;
		break;
	case 3:
		s = 1./6.;
		break;
	default:
        throw MechanicsException ( "[NuTo::RungeKutta4::GetStageWeights] rStage>3 not implemented." );
	}
	return s;
}
#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::RungeKutta4::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::RungeKutta4::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::RungeKutta4::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::RungeKutta4::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::RungeKutta4::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::RungeKutta4::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::RungeKutta4::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of RungeKutta4" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RungeKuttaBase);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of RungeKutta4" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::RungeKutta4::GetTypeId()const
{
    return std::string("RungeKutta4");
}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::RungeKutta4::Restore (const std::string &filename, std::string rType )
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
                throw MechanicsException ( "[RungeKutta4::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[RungeKutta4::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[RungeKutta4::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
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
        throw MechanicsException ( "[RungeKutta4::Restore]Unhandled exception." );
    }
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::RungeKutta4::Save (const std::string &filename, std::string rType )const
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
            throw MechanicsException ( "[RungeKutta4::Save]File type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception& e )
    {
        std::string s ( std::string ( "[RungeKutta4::Save]File save exception in boost - " ) +std::string ( e.what() ) );
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
        throw MechanicsException ( "[RungeKutta4::Save] Unhandled exception." );
    }
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::RungeKutta4)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
