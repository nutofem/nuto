// $Id: Newmark.cpp 575 2011-09-20 18:05:35Z unger3 $

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

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverMKLPardiso.h"

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/NewmarkIndirect.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/mechanics/constraints/ConstraintLinearEquation.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::NewmarkIndirect::NewmarkIndirect (StructureBase& rStructure)  : NewmarkBase (rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::NewmarkIndirect::Info()const
{
	NewmarkBase::Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NewmarkIndirect::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NewmarkIndirect::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NewmarkIndirect::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NewmarkIndirect::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NewmarkIndirect::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NewmarkIndirect::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NewmarkIndirect::serialize(Archive & ar, const unsigned int version)
{
	#ifdef DEBUG_SERIALIZATION
	    std::cout << "start serialization of Newmark" << "\n";
	#endif
	    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NewmarkBase);


    #ifdef DEBUG_SERIALIZATION
         std::cout << "finish serialization of Newmark" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::Error::eError NuTo::NewmarkIndirect::Solve(StructureBase& rStructure, double rTimeDelta)
{
    throw MechanicsException("[NuTo::NewmarkIndirect::Solve] to be implemented.");
	return NuTo::Error::SUCCESSFUL;
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::NewmarkIndirect::GetTypeId()const
{
    return std::string("NewmarkIndirect");
}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::NewmarkIndirect::Restore (const std::string &filename, std::string rType )
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
				throw MechanicsException ( "[Newmark::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
		}
		else if (rType=="XML")
		{
			boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
			oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
			if ( tmpString!=GetTypeId() )
				throw MechanicsException ( "[Newmark::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
		}
		else if (rType=="TEXT")
		{
			boost::archive::text_iarchive ota ( ifs, std::ios::binary );
			ota & boost::serialization::make_nvp ( "Object_type", tmpString );
			if ( tmpString!=GetTypeId() )
				throw MechanicsException ( "[Newmark::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            ota & boost::serialization::make_nvp(tmpString.c_str(), *this);
		}
		else
		{
			throw MathException ( "[NewmarkIndirect::Restore]File type not implemented" );
		}
	}
	catch ( MechanicsException &e )
	{
		throw e;
	}
	catch ( std::exception &e )
	{
		throw MechanicsException ( e.what() );
	}
	catch ( ... )
	{
		throw MechanicsException ( "[NewmarkIndirect::Restore]Unhandled exception." );
	}
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::NewmarkIndirect::Save (const std::string &filename, std::string rType )const
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
			throw MechanicsException ( "[NewmarkIndirect::Save]File type not implemented." );
		}
	}
	catch ( boost::archive::archive_exception e )
	{
		std::string s ( std::string ( "[NewmarkIndirect::Save]File save exception in boost - " ) +std::string ( e.what() ) );
		std::cout << s << "\n";
		throw MathException ( s );
	}
	catch ( MechanicsException &e )
	{
		throw e;
	}
	catch ( std::exception &e )
	{
		throw MechanicsException ( e.what() );
	}
	catch ( ... )
	{
		throw MechanicsException ( "[NewmarkIndirect::Save]Unhandled exception." );
	}
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NewmarkIndirect)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
