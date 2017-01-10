// $Id: NystroemQinZhu.cpp 575 2011-09-20 18:05:35Z unger3 $

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

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/NystroemQinZhu.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"


//! @brief constructor
//! @param mDimension number of nodes
NuTo::NystroemQinZhu::NystroemQinZhu (StructureBase* rStructure)  : NystroemBase (rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::NystroemQinZhu::Info()const
{
	TimeIntegrationBase::Info();
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous step)
// so essentially it's c_n-c_(n-1)
double NuTo::NystroemQinZhu::GetStageTimeFactor(int rStage)const
{
	assert(rStage<3);
	double s;
	switch(rStage)
	{
	case 0:
		s = (3.+sqrt(3.))/6.;
		break;
	case 1:
		s = (3.-sqrt(3.))/6.;
		break;
	case 2:
		s = (3.+sqrt(3.))/6.;
		break;
	default:
		throw MechanicsException("[NuTo::NystroemQinZhu::GetStageTimeFactor] error with stage number.");
	}
	return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous step)
// so essentially it's c_n-c_(n-1)
bool NuTo::NystroemQinZhu::HasTimeChanged(int rStage)const
{
	assert(rStage<3);
	bool s;
	switch(rStage)
	{
	case 0:
		s = true;
		break;
	case 1:
		s = true;
		break;
	case 2:
		s = true;
		break;
	default:
		throw MechanicsException("[NuTo::NystroemQinZhu::HasTimeChanged] error with stage number.");
	}
	return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::NystroemQinZhu::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage)const
{
	assert(rStage<3);
	assert(rWeight.size()==2);
	switch(rStage)
	{
	case 0:
		break;
	case 1:
		rWeight[0] = (2.-sqrt(3.))/12.;
		break;
	case 2:
		rWeight[0] = 0.0;
		rWeight[1] = sqrt(3.)/6.;
		break;
	default:
		throw MechanicsException("[NuTo::NystroemQinZhu::GetStageDerivativeFactor] error with stage number.");
	}
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::NystroemQinZhu::GetStageWeights1(int rStage)const
{
	assert(rStage<3);
	double s;
	switch(rStage)
	{
	case 0:
		s = (5.-3.*sqrt(3.))/24.;
		break;
	case 1:
		s = (3.+sqrt(3.))/12.;
		break;
	case 2:
		s = (1.+sqrt(3.))/24.;
		break;
	default:
		throw MechanicsException("[NuTo::NystroemQinZhu::GetStageWeights1] error with stage number.");
	}
	return s;
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::NystroemQinZhu::GetStageWeights2(int rStage)const
{
	assert(rStage<3);
	double s;
	switch(rStage)
	{
	case 0:
		s = (3.-2.*sqrt(3.))/12.;
		break;
	case 1:
		s = 0.5;
		break;
	case 2:
		s = (3.+2.*sqrt(3.))/12.;
		break;
	default:
		throw MechanicsException("[NuTo::NystroemQinZhu::GetStageWeights2] error with stage number.");
	}
	return s;
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
//! check paper "Lax-Wendroff and Nystrom methods for seismic modelling" by Jing-Bo Chen
double NuTo::NystroemQinZhu::CalculateCriticalTimeStep()const
{
	double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
	return 2.58651889452/std::sqrt(maxGlobalEigenValue);
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NystroemQinZhu::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NystroemQinZhu::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NystroemQinZhu::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NystroemQinZhu::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NystroemQinZhu::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NystroemQinZhu::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NystroemQinZhu::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of NystroemQinZhu" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NystroemBase);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of NystroemQinZhu" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::NystroemQinZhu::GetTypeId()const
{
    return std::string("NystroemQinZhu");
}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::NystroemQinZhu::Restore (const std::string &filename, std::string rType )
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
                throw MechanicsException ( "[NystroemQinZhu::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[NystroemQinZhu::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[NystroemQinZhu::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
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
        throw MechanicsException ( "[NystroemQinZhu::Restore]Unhandled exception." );
    }
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::NystroemQinZhu::Save (const std::string &filename, std::string rType )const
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
            throw MechanicsException ( "[NystroemQinZhu::Save]File type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception& e )
    {
        std::string s ( std::string ( "[NystroemQinZhu::Save]File save exception in boost - " ) +std::string ( e.what() ) );
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
        throw MechanicsException ( "[NystroemQinZhu::Save] Unhandled exception." );
    }
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NystroemQinZhu)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
