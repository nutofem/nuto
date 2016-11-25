// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION


#include "nuto/mechanics/timeIntegration/NewmarkBase.h"


namespace NuTo
{

class StructureOutputBlockMatrix;
class StructureOutputBlockVector;

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

    void SetPerformLineSearch(bool rPerformLineSearch)
    {
        mPerformLineSearch = rPerformLineSearch;
    }

    bool GetPerformLineSearch() const
    {
        return mPerformLineSearch;
    }

    int GetVerboseLevel() const
    {
        return mVerboseLevel;
    }

    void SetVerboseLevel(int rVerboseLevel)
    {
        mVerboseLevel = rVerboseLevel;
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
#endif // ENABLE_SERIALIZATION

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    virtual NuTo::eError Solve(double rTimeDelta);

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    virtual std::string GetTypeId()const;

protected:
    StructureOutputBlockVector CalculateDof1(
            const StructureOutputBlockVector& rDeltaDof_dt0,
            const StructureOutputBlockVector& rDof_dt1,
            const StructureOutputBlockVector& rDof_dt2, double rTimeStep) const;

    StructureOutputBlockVector CalculateDof2(
            const StructureOutputBlockVector& rDeltaDof_dt0,
            const StructureOutputBlockVector& rDof_dt1,
            const StructureOutputBlockVector& rDof_dt2, double rTimeStep) const;

    //! @brief ... builds the modified hessian matrix (including cmat) and solves the system
    //! @param rHessian_dt0 ... hessian_dt0 matrix
    //! @param rHessian_dt1 ... hessian_dt1 matrix
    //! @param rHessian_dt2 ... hessian_dt2 matrix
    //! @param rResidualMod ... modified residual (including the cmatrix)
    //! @param rTimeStep ... current time step
    //! @return ... deltaDof_dt0.J
    //! @remark ... If hessian_dt0 is constant, its values are preserved (hessianMod = temporary matrix). Otherwise, hessian0 will be used (hessianMod = hessian0)
    BlockFullVector<double> BuildHessianModAndSolveSystem(
                  StructureOutputBlockMatrix& rHessian_dt0,
            const StructureOutputBlockMatrix& rHessian_dt1,
            const StructureOutputBlockMatrix& rHessian_dt2,
            const BlockFullVector<double>& rResidualMod, double rTimeStep) const;

    StructureOutputBlockVector CalculateResidual(
            const StructureOutputBlockVector& rIntForce,
            const StructureOutputBlockVector& rExtForce,
            const StructureOutputBlockMatrix& rHessian2,
            const StructureOutputBlockVector& rDof_dt1,
            const StructureOutputBlockVector& rDof_dt2) const;

    //! @brief Calculates (if needed) the residual.K part for the post-processing. Since it is not needed for the actual time integration
    //! its calculation is skipped if Cmat has only zero entries.
    void CalculateResidualKForPostprocessing(
            StructureOutputBlockVector& rResidual,
            const StructureOutputBlockMatrix& rHessian_dt2,
            const StructureOutputBlockVector& rDof_dt1,
            const StructureOutputBlockVector& rDof_dt2) const;

    void CalculateMuDampingMatrix(StructureOutputBlockMatrix& rHessian_dt1, const StructureOutputBlockMatrix& rHessian_dt2) const;

    void CalculateResidualTrial(
            StructureOutputBlockVector& rResidual,
            const BlockFullVector<double>& rDeltaBRHS,
            const StructureOutputBlockMatrix& rHessian_dt0,
            const StructureOutputBlockMatrix& rHessian_dt1,
            const StructureOutputBlockMatrix& rHessian_dt2,
            const StructureOutputBlockVector& rDof_dt1,
            const StructureOutputBlockVector& rDof_dt2,
            double rTimeStep) const;


    //! @brief Prints Info about the current calculation stage
    void PrintInfoStagger() const;

    //! @brief Prints Info about the current iteration
    void PrintInfoIteration(const BlockScalar &rNormResidual, int rIteration) const;

protected:

    double mMinLineSearchStep;

    bool mPerformLineSearch;

    int mVisualizeResidualTimeStep;

    //! @brief Controls the output verbosity (0 = silent)
    int mVerboseLevel = 1;


};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::NewmarkDirect)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



