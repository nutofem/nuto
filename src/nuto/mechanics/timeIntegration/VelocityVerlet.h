// $Id$

#ifndef VelocityVerlet_H
#define VelocityVerlet_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (Newmark, but you can use it for statics as well with setting the flag isDynamic to false)
class VelocityVerlet : public TimeIntegrationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    VelocityVerlet(StructureBase* rStructure);

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif// SWIG

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    //! @brief ... save the object to a file
    void Restore (const std::string &filename, std::string rType );

	//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
	void Save (const std::string &filename, std::string rType )const;
#endif // ENABLE_SERIALIZATION

    //! @brief perform the time integration
    //! @param rStructure ... structure
    //! @param rTimeDelta ... length of the simulation
    NuTo::Error::eError Solve(double rTimeDelta);

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const;

    //! @brief ... Return time step
    double GetTimeStep()const
    {
    	return mTimeStep;
    }

    //! @brief ... Return time step
    void SetTimeStep(double rTimeStep)
    {
    	mTimeStep = rTimeStep;
    }

protected:
    //empty private construct required for serialization
    VelocityVerlet(){};
	//time step for the time integration, be careful not to make it smaller than the critical time step
    double mTimeStep;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::VelocityVerlet)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



#endif // VelocityVerlet_H
