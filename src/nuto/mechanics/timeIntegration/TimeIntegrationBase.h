// $Id: TimeIntegrationBase.h 593 2012-02-03 23:50:23Z unger3 $

#ifndef TIMEINTEGRATIONBASE_H
#define TIMEINTEGRATIONBASE_H

#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/base/NuToObject.h"
#include "nuto/base/ErrorEnum.h"

namespace NuTo
{
class StructureBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all mechanical structures
class TimeIntegrationBase : public NuToObject
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    TimeIntegrationBase();

    //! @brief deconstructor
    virtual ~TimeIntegrationBase()
    {}

    //! @brief deconstructor
    virtual NuTo::Error::eError Solve(NuTo::StructureBase* rStructure)const=0;

    //! @brief sets the constraint equation whose RHS is incrementally increased in each load step / time step
    void SetConstraintLoad(int rConstraintLoad)
    {
    	mConstraintLoad = rConstraintLoad;
    }

    //! @brief sets the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
    void SetConstraintRHSDelta(int rConstraintRHSDelta)
    {
    	mConstraintRHSDelta = rConstraintRHSDelta;
    }

    //! @brief sets the delta time
    void SetDeltaTime(int rTimeDelta)
    {
    	mTimeDelta = rTimeDelta;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info()const;
protected:
    int mConstraintLoad;
	double mConstraintRHSDelta;
	double mTime;
	double mTimeDelta;
    int mLoadStep;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::TimeIntegrationBase)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // TIMEINTEGRATIONBASE_H
