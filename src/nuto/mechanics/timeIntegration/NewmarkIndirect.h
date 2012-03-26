// $Id$

#ifndef NEWMARKINDIRECT_H
#define NEWMARKINDIRECT_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/NewmarkBase.h"

namespace NuTo
{
class ConstraintLinearEquation;
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (static Newton Raphson or Newmark for dynamics)
class NewmarkIndirect : public NewmarkBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NewmarkIndirect();

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
    virtual void Restore (const std::string &filename, std::string rType );

	//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
	virtual void Save (const std::string &filename, std::string rType )const;
#endif // ENABLE_SERIALIZATION

    //! @brief perform the time integration
    //! @param rStructure ... structure
    //! @param rTimeDelta ... length of the simulation
    NuTo::Error::eError Solve(StructureBase& rStructure, double rTimeDelta);

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    virtual std::string GetTypeId()const;

protected:
    ConstraintLinearEquation* mLinearConstraint;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::NewmarkIndirect)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



#endif // NEWMARKINDIRECT_H
