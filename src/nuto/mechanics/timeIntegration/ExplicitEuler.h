// $Id$

#ifndef ExplicitEuler_H
#define ExplicitEuler_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"

namespace NuTo
{
class ExplicitEuler : public TimeIntegrationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    ExplicitEuler(StructureBase* rStructure);

    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    bool HasCriticalTimeStep()const
    {
    	return true;
    }

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    double CalculateCriticalTimeStep()const;

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif// SWIG

    //! @brief ... restore the object from a ExplicitEulere
    //! @param ExplicitEulerename ... ExplicitEulerename
    //! @param aType ... type of ExplicitEulere, either BINARY, XML or TEXT
    //! @brief ... save the object to a ExplicitEulere
    void Restore (const std::string &ExplicitEulerename, std::string rType );

	//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
    //! @param ExplicitEulerename ... ExplicitEulerename
    //! @param aType ... type of ExplicitEulere, either BINARY, XML or TEXT
    void Save (const std::string &ExplicitEulerename, std::string rType )const;
#endif // ENABLE_SERIALIZATION

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the ExplicitEulere
    //!            in case of restoring from a ExplicitEulere with the wrong object type, the ExplicitEulere id is printed
    //! @return    class name
    std::string GetTypeId()const;

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    NuTo::Error::eError Solve(double rTimeDelta);

protected:
    //empty private construct required for serialization
    ExplicitEuler(){}
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::ExplicitEuler)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



#endif
