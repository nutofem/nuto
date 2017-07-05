#ifdef _OPENMP
#include <omp.h>
#endif

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/RungeKutta2.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"


//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKutta2::RungeKutta2(StructureBase* rStructure)
    : RungeKuttaBase(rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKutta2::Info() const
{
    TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::RungeKutta2::CalculateCriticalTimeStep() const
{
    double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.0 / std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
double NuTo::RungeKutta2::GetStageTimeFactor(int rStage) const
{
    assert(rStage < 2);
    double s;
    switch (rStage)
    {
    case 0:
        s = 0.;
        break;
    case 1:
        s = 0.5;
		break;
	default:
        throw Exception ( "[NuTo::RungeKutta2::GetStageTimeFactor] rStage>3 not implemented." );
	}
	return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous
//! step)
// so essentially it's c_n-c_(n-1)
bool NuTo::RungeKutta2::HasTimeChanged(int rStage) const
{
    assert(rStage < 2);
    bool s;
    switch (rStage)
    {
    case 0:
        s = false; // same as last step from the last iteration
        break;
    case 1:
		s = true;
		break;
	default:
        throw Exception ( "[NuTo::RungeKutta2::HasTimeChanged] rStage>3 not implemented." );
	}
	return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::RungeKutta2::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const
{
    assert(rStage < 2);
    assert(rWeight.size() == 1);
    switch (rStage)
    {
    case 0:
        break;
    case 1:
        rWeight[0] = 0.5;
		break;
	default:
        throw Exception ( "[NuTo::RungeKutta2::GetStageDerivativeFactor] rStage>3 not implemented." );
	}
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::RungeKutta2::GetStageWeights(int rStage) const
{
    assert(rStage < 2);
    double s;
    switch (rStage)
    {
    case 0:
        s = 0.0;
        break;
    case 1:
        s = 1.0;
		break;
	default:
        throw Exception ( "[NuTo::RungeKutta2::GetStageWeights] rStage>3 not implemented." );
	}
	return s;
}
#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::RungeKutta2::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::RungeKutta2::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::RungeKutta2::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::RungeKutta2::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::RungeKutta2::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::RungeKutta2::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::RungeKutta2::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of RungeKutta2" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RungeKuttaBase);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of RungeKutta2" << "\n";
    #endif
}


//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::RungeKutta2::Restore (const std::string &filename, std::string rType )
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
                throw Exception ( "[RungeKutta2::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw Exception ( "[RungeKutta2::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw Exception ( "[RungeKutta2::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            ota & boost::serialization::make_nvp(tmpString.c_str(), *this);
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
        throw Exception ( "[RungeKutta2::Restore]Unhandled exception." );
    }
    return s;
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::RungeKutta2::Save (const std::string &filename, std::string rType )const
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        std::string tmpStr ( GetTypeId() );
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
            throw Exception ( "[RungeKutta2::Save]File type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception& e )
    {
        std::string s ( std::string ( "[RungeKutta2::Save]File save exception in boost - " ) +std::string ( e.what() ) );
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
        throw Exception ( "[RungeKutta2::Save] Unhandled exception." );
    }
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::RungeKutta2)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
