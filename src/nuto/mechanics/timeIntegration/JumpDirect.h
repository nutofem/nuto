// $Id$

#ifndef JUMPDIRECT_H
#define JUMPDIRECT_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/NewmarkBase.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include <boost/math/constants/constants.hpp>
#include <math.h>

namespace NuTo
{
//! @author Vitaliy M. Kindrachuk
//! @date December 2014
//! @brief ... standard class for implicit time integration including cycle jump based on a Fourier formulation
class JumpDirect : public NewmarkDirect
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    JumpDirect(StructureBase* rStructure);

    void SetMinLineSearchStep(double rMinLineSearchStep)
    {
    	mMinLineSearchStep = rMinLineSearchStep;
    }

    double GetMinLineSearchStep()const
    {
    	return mMinLineSearchStep;
    }

    //!@brief sets the number of increments of a single cycle, this should be divisible by 4, because the cycle consist of 4 parts
    //!@param sets the number of increments of a single cycle, this should be divisible by 4, because the cycle consist of 4 parts
    void SetHarmonicIncrementation(int rHarmonicIncrementation)
    {
    	if (fmod(rHarmonicIncrementation, 4) == 0) {
    		mHarmonicIncrementation = rHarmonicIncrementation;
    	} else {
    		mHarmonicIncrementation = (rHarmonicIncrementation - fmod(rHarmonicIncrementation,4) + 4);
    	}
    }

    //!@brief gets the number of increments of a single cycle, this should be divisible by 4, because the cycle consist of 4 parts
    //!@param gets the number of increments of a single cycle, this should be divisible by 4, because the cycle consist of 4 parts
    int GetHarmonicIncrementation()const
    {
    	return mHarmonicIncrementation;
    }

    //! @brief if mHarmonicExtrapolation == false then performs a Newmark, otherwise performs Newmark till 3 cycles and then cycle jump
    //! @param if mHarmonicExtrapolation == false then performs a Newmark, otherwise performs Newmark till 3 cycles and then cycle jump
    void SetHarmonicExtrapolation(bool rHarmonicExtrapolation)
    {
    	if (mHarmonicExcitation == true) {
    		mHarmonicExtrapolation = rHarmonicExtrapolation;
		} else {
			throw MechanicsException("[NuTo::JumpDirect::SetHarmonicConstraint] SetHarmonicConstraint at first!");
		}

    }

    //! @brief if mHarmonicExtrapolation == false then performs a Newmark, otherwise performs Newmark till 3 cycles and then cycle jump
    //! @param if mHarmonicExtrapolation == false then performs a Newmark, otherwise performs Newmark till 3 cycles and then cycle jump
    bool GetHarmonicExtrapolation()const
    {
    	return mHarmonicExtrapolation;
    }

    //!@brief sets the harmonic excitations in the format {DisplacementAmplitude,Frequency,NumberOfCycles}
    //!@param sets the harmonic excitations in the format {DisplacementAmplitude,Frequency,NumberOfCycles}
    // rHarmonicConstraintFactor[*,0] = displacement amplitude, arbitrary sign
    // rHarmonicConstraintFactor[*,1] = frequency, positive
    // rHarmonicConstraintFactor[*,2] = number of cycles, at least 1
    // The harmonic excitations are applied to the mTimeDependentConstraint !
    // Currently is only a monoharmonic constraint applied, that is rHarmonicConstraintFactor.GetNumRows()=1
    void SetHarmonicConstraint(const NuTo::FullMatrix<double,Eigen::Dynamic,3>& rHarmonicConstraintFactor)
    {
    	if (mTimeDependentConstraint == -1) {
    		throw MechanicsException("[NuTo::JumpDirect::SetHarmonicConstraint] the harmonic excitation can be only applied to the time dependent constraint.");
		}

    	if (rHarmonicConstraintFactor.GetNumRows()>1)
    		throw MechanicsException("[NuTo::JumpDirect::SetHarmonicConstraint] currently a monoharmonic constraint is implemented, number of rows must be 1.");

    	//check, whether frequencies and number of cycles are positive
    	for (int count=0; count<rHarmonicConstraintFactor.GetNumRows(); count++)
    	{
    		if (rHarmonicConstraintFactor(count,1)<=0.)
    			throw MechanicsException("[NuTo::JumpDirect::SetHarmonicConstraint] the frequency should always be positive.");
    		if (rHarmonicConstraintFactor(count,2)<=1.)
    			throw MechanicsException("[NuTo::JumpDirect::SetHarmonicConstraint] number of cycles should be at least 1.");
    	}

    	mHarmonicExcitation = true;
    	mHarmonicConstraintFactor = rHarmonicConstraintFactor;
    }

    //! @brief returns true, if there is a harmonic excitation applied
    bool GetHarmonicExcitation()const
    {
    	return mHarmonicExcitation;
    }

    //! @brief apply the new rhs of the constraints as a function of the current time delta
    //! @brief the time dependent constraint is followed by a sine excitation
    double CalculateTimeDependentConstraintFactor(double curTime)
    {
    	if (mTimeDependentConstraintFactor.GetNumRows()==0)
    		throw MechanicsException("[NuTo::JumpDirect::CalculateTimeDependentConstraintFactor] the harmonic excitation can be only applied to the time dependent constraint.");

    	// calculate the end time of the time dependent constraint. Since this time the constraint is harmonic
    	double timeDependentConstraintTime(mTimeDependentConstraintFactor(mTimeDependentConstraintFactor.GetNumRows()-1,0));
		double timeDependentConstraintFactor;

    	if (curTime <= timeDependentConstraintTime) {
        	// calculate timeDependentConstraintFactor during the time dependent constraint
    		return NewmarkDirect::CalculateTimeDependentConstraintFactor(curTime);
		} else {
			// time dependent constraint is continued by the harmonic excitation

			// calculate constraint at the end time of the time dependent constraint
			timeDependentConstraintFactor = NewmarkDirect::CalculateTimeDependentConstraintFactor(timeDependentConstraintTime);

			if (mHarmonicConstraintFactor.GetNumRows()==0) {
				return timeDependentConstraintFactor;
			}

			// add harmonic excitation
			double frequency(mHarmonicConstraintFactor(0,1));
			double amplitude(mHarmonicConstraintFactor(0,0));
			const double pi = boost::math::constants::pi<double>();
			timeDependentConstraintFactor += amplitude*sin(2*pi*frequency*(curTime - timeDependentConstraintTime));
			return timeDependentConstraintFactor;
		}
    }

    void SetTimeAndTimeStep(double &curTime, double &timeStep, double rTimeDelta)
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

        // calculate time step during harmonic excitation, the number of increments should be divisible by 4
        if (mHarmonicExcitation == true && fabs(timeStep - 1./(mHarmonicIncrementation*mHarmonicConstraintFactor(0,1))) > 0.001*timeStep) {
        	curTime += 1./(mHarmonicIncrementation*mHarmonicConstraintFactor(0,1))-timeStep;
        	timeStep = 1./(mHarmonicIncrementation*mHarmonicConstraintFactor(0,1));
        	return;
		}
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

    //! @brief performs the time integration
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
    JumpDirect(){};
	double mMinLineSearchStep;
	bool mHarmonicExcitation;
	bool mHarmonicExtrapolation;
	int mHarmonicIncrementation;
	NuTo::FullMatrix<double,Eigen::Dynamic,3> mHarmonicConstraintFactor;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::JumpDirect)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



#endif // JumpDirect_H
