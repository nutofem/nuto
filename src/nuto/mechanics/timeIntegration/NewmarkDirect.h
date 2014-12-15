// $Id$

#ifndef NEWMARKDIRECT_H
#define NEWMARKDIRECT_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/NewmarkBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (static Newton Raphson or NewmarkDirect for dynamics)
class NewmarkDirect : public NewmarkBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NewmarkDirect(StructureBase* rStructure);

    void SetMinLineSearchStep(double rMinLineSearchStep)
    {
    	mMinLineSearchStep = rMinLineSearchStep;
    }

    double GetMinLineSearchStep()const
    {
    	return mMinLineSearchStep;
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

    virtual void SetTimeAndTimeStep(double &curTime, double &timeStep, double rTimeDelta)
    {
    	// calculate time and time step close to the end time of the integration interval
        if (rTimeDelta-curTime<0.5*timeStep)
        {
            timeStep += rTimeDelta-curTime;
            curTime = rTimeDelta;
            return;
        }

    	// calculate time and time step during the time dependent constraint (where the linear interpolation is performed)
    	if (mTimeDependentConstraintFactor.GetNumRows()!=0 && mAutomaticTimeStepping)
    	{
    		int curStep(0);
    		while (mTimeDependentConstraintFactor(curStep,0)<curTime && curStep<mTimeDependentConstraintFactor.GetNumRows()-1)
    			curStep++;
    		if (curStep==0)
    			curStep++;

    		//extract the two data points
    		double t1 = mTimeDependentConstraintFactor(curStep-1,0);
			double t2 = mTimeDependentConstraintFactor(curStep,0);

			// if curTime is close to the end of the time increment then set it to the end
            if (t2-curTime<0.5*timeStep && curTime<t2)
            {
                timeStep += t2-curTime;
                curTime = t2;
                return;
            }
            // if curTime jumps over a data point, then set it to this data point
            if (curTime-timeStep < t1 - 0.2*timeStep) {
    			timeStep -= curTime -t1;
    			curTime = t1;
    			return;
    		}
            // if curTime jumps over the end of the time dependent constraint, then set it
            // to this end (this point is at the same time the beginning of the harmonic excitation)
            if (curTime-timeStep < t2 - 0.2*timeStep && curTime > t2) {
    			timeStep -= curTime -t2;
    			curTime = t2;
    			return;
    		}
            if (curTime<t2) {
				return;
			}
		}

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
    NewmarkDirect(){};
	double mMinLineSearchStep;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::NewmarkDirect)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



#endif // NewmarkDirect_H
