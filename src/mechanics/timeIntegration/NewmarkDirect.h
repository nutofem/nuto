// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION


#include "mechanics/timeIntegration/NewmarkBase.h"


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
    bool HasCriticalTimeStep()const override
    {
    	return false;
    }

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    double CalculateCriticalTimeStep()const  override
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
    virtual void Solve(double rTimeDelta) override;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const override;

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



