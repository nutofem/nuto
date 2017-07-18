#pragma once

#include <vector>

#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/dofSubMatrixSolvers/SolverBase.h"
#include "mechanics/timeIntegration/TimeControl.h"

namespace NuTo
{
class CallbackInterface;
class NodeBase;
class StructureBase;
class TimeDependencyBase;
class PostProcessor;
enum class eError;

template <typename T>
class BlockFullVector;

namespace IpData
{
enum class eIpStaticDataType;
} // namespace IpData

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all mechanical structures
class TimeIntegrationBase
{
public:
    //! @brief constructor
    TimeIntegrationBase(StructureBase* rStructure);

    //! @brief deconstructor
    virtual ~TimeIntegrationBase();

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    virtual void Solve(double rTimeDelta) = 0;

    //! @brief sets the delta rhs of the constrain equation whose RHS is incrementally increased in each load step /
    //! time step
    void ResetForNextLoad();

    //! @brief Gets the Blockscalar with the set residual tolerances
    //! @return Blockscalar with the set residual tolerances
    const BlockScalar& GetToleranceResidual() const;

    //! @brief Sets the residual tolerance for a specific DOF
    //! param rDof: degree of freedom
    //! param rTolerance: tolerance
    void SetToleranceResidual(Node::eDof rDof, double rTolerance);

    //! @brief Updates the Rhs for all constraints
    //! @param rCurrentTime ... current time
    //! @remark remove the second argument rDof
    void UpdateConstraints(double rCurrentTime);

    //! @brief sets a scalar time dependent multiplication factor for the external loads
    //! @param rLoadRHSFactor ... first row time, second row scalar factor to calculate the external load (linear
    //! interpolation in between,  afterwards linear extrapolation)
    void SetTimeDependentLoadCase(int rLoadCase, const Eigen::MatrixXd& rLoadRHSFactor);

    //! @brief calculate the external force vector (mStatic and m) as a function of time delta
    virtual void CalculateStaticAndTimeDependentExternalLoad();

    //! @brief calculate the current external force as a function of time delta
    //! @param curTime ... current time in the load step
    //! @return ... external load vector
    virtual StructureOutputBlockVector CalculateCurrentExternalLoad(double curTime);

    //! @brief calculates the norm of the residual, can include weighting
    //! @param rResidual ... residual
    double CalculateNorm(const BlockFullVector<double>& rResidual) const;

    //! @brief sets the  time step for the time integration procedure (initial value)
    virtual void SetTimeStep(double rTimeStep)
    {
        mTimeControl.SetTimeStep(rTimeStep);
    }

    //! @brief returns the  time step for the time integration procedure (current value)
    virtual double GetTimeStep() const
    {
        return mTimeControl.GetTimeStep();
    }
    //! @brief sets the maximum time step for the time integration procedure
    void SetMaxTimeStep(double rMaxTimeStep)
    {
        mTimeControl.SetMaxTimeStep(rMaxTimeStep);
    }


    //! @brief returns the maximum time step for the time integration procedure
    double GetMaxTimeStep() const
    {
        return mTimeControl.GetMaxTimeStep();
    }

    //! @brief sets the minimum time step for the time integration procedure
    void SetMinTimeStep(double rMinTimeStep)
    {
        mTimeControl.SetMinTimeStep(rMinTimeStep);
    }

    //! @brief returns the g
    std::vector<int> GetVecGroupNodesReactionForces() const
    {
        return mVecGroupNodesReactionForces;
    }

    //! @brief getter for iteration count, initialized to zero, incremented on each iteration.
    int GetNumIterations() const
    {
        return mIterationCount;
    }

    //! @brief sets automatic time stepping (on or off)
    void SetAutomaticTimeStepping(bool rAutomaticTimeStepping)
    {
        mAutomaticTimeStepping = rAutomaticTimeStepping;
        if (mAutomaticTimeStepping)
        {
            mTimeControl.UseDefaultAutomaticTimestepping();
        }
        else
        {
            mTimeControl.UseEquidistantTimestepping();
        }
    }


    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    virtual bool HasCriticalTimeStep() const = 0;

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    virtual double CalculateCriticalTimeStep() const = 0;

    void ConnectCallback(CallbackInterface* rCallback)
    {
        mCallback = rCallback;
    }


    //! @brief ... Adds a calculation step to each timestep
    //! param rActiveDofs ... active Dofs of the calculation step
    void AddCalculationStep(const std::set<Node::eDof>& rActiveDofs);

    //! @brief ... Sets the number of calculation steps per timestep
    //! param rNumSteps ... number of calculation steps per timestep
    void SetNumCalculationSteps(int rNumSteps);

    //! @brief ... Sets the active Dofs of a calculation step
    //! param rStepNum ... step number
    //! param rActiveDofs ... active Dofs of the calculation step
    void SetActiveDofsCalculationStep(int rStepNum, const std::set<Node::eDof>& rActiveDofs);

#ifndef SWIG
    //! @brief Sets the solver
    void SetSolver(std::unique_ptr<SolverBase> rSolver)
    {
        mSolver = std::move(rSolver);
    }
#endif

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info() const;

    bool GetShowTime() const;

    void SetShowTime(bool showTime);

    const PostProcessor& PostProcessing() const
    {
        return *mPostProcessor;
    }

    PostProcessor& PostProcessing()
    {
        return *mPostProcessor;
    }

protected:
    //! @brief extracts all dof values
    //! @return ret[0] are the DOF values, ret[1], ret[2] the 1st and 2nd time derivative, respectively
    std::array<StructureOutputBlockVector, 3> ExtractDofValues() const;

    const BlockFullVector<double>& UpdateAndGetConstraintRHS(double rCurrentTime);

    const BlockFullVector<double>& UpdateAndGetAndMergeConstraintRHS(double rCurrentTime,
                                                                     StructureOutputBlockVector& rDof_dt0);

    StructureBase* mStructure; //!< structure belonging to the time integration scheme

    std::unique_ptr<SolverBase> mSolver; //!< sparse matrix solver

    int mTimeDependentLoadCase = -1; //!< time dependent load case number

    // includes for each time step the scalar factor for the load case
    // the time step is given relative to mTimeDelta
    Eigen::MatrixXd mTimeDependentLoadFactor;

    NuTo::StructureOutputBlockVector mLoadVectorStatic; //!< static external load vector
    NuTo::StructureOutputBlockVector mLoadVectorTimeDependent; //!< dynamic external load vector

    bool mAutomaticTimeStepping = false; //!< adapt the time step based on the number of iterations required (or
    //! decrease, if no convergence can be achieved)

    bool mCheckCoefficientMatrix = false; //!< if set to true, checks the coefficient matrix in each sub step

    BlockScalar mToleranceResidual;

    int mIterationCount = 0; //!< iteration count


    int mLoadStep = 0; //!< load step number is increased after each converged step (used for successive output)

    std::vector<int> mVecGroupNodesReactionForces; //!< vector of groups of nodes for which the residual (corresponding
    //! to the reaction forces induced by constraints) is given as output

    CallbackInterface* mCallback;

    std::vector<std::set<Node::eDof>> mStepActiveDofs; //!< stores wich Dofs are active in which calculation step

    bool mShowTime = true;
    TimeControl mTimeControl;

    std::unique_ptr<PostProcessor> mPostProcessor;
};
} // namespace NuTo
