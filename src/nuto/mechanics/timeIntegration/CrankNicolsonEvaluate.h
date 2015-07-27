// $Id$

#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, BAM
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (static Newton Raphson or CrankNicolsonEvaluate for dynamics)
class CrankNicolsonEvaluate : public TimeIntegrationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    CrankNicolsonEvaluate(StructureBase* rStructure);

    void SetToleranceForce(double rToleranceForce)
    {
    	mToleranceForce = rToleranceForce;
    }

    double GetToleranceForce()const
    {
    	return mToleranceForce;
    }

    void SetMinLineSearchStep(double rMinLineSearchStep)
    {
    	mMinLineSearchStep = rMinLineSearchStep;
    }

    double GetMinLineSearchStep()const
    {
    	return mMinLineSearchStep;
    }

    void SetPerformLineSearch(bool rPerformLineSearch)
    {
        mPerformLineSearch = rPerformLineSearch;
    }

    bool GetPerformLineSearch() const
    {
        return mPerformLineSearch;
    }

    void SetMaxNumIterations(unsigned int rMaxNumIterations)
    {
        mMaxNumIterations = rMaxNumIterations;
    }

    void SetCheckEquilibriumOnStart(bool rCheckEquilibriumOnStart)
    {
        mCheckEquilibriumOnStart = rCheckEquilibriumOnStart;
    }

    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    bool HasCriticalTimeStep()const
    {
    	return false;
    }

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    double CalculateCriticalTimeStep()const
    {
    	return 0;
    }

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
    //! @param rTimeDelta ... length of the simulation
    NuTo::Error::eError Solve(double rTimeDelta);

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    virtual std::string GetTypeId()const;

protected:
    //empty private construct required for serialization
    CrankNicolsonEvaluate(){};
	double mMinLineSearchStep;

	bool mPerformLineSearch;

    bool mCheckEquilibriumOnStart = true;

	int mVisualizeResidualTimeStep;

	double mToleranceForce;
	int mMaxNumIterations;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::CrankNicolsonEvaluate)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



#endif // CRANK_NICOLSON_H
