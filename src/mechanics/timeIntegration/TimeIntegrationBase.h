#pragma once

#include <vector>
#include <boost/ptr_container/ptr_map.hpp>

// member
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/dofSubMatrixSolvers/SolverBase.h"

namespace NuTo
{
class CallbackInterface;
class NodeBase;
class ResultBase;
class StructureBase;
class TimeDependencyBase;
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

    //! @brief extracts all dof values
    //! @param rDof_dt0 ... 0th time derivative
    //! @param rDof_dt1 ... 1st time derivative
    //! @param rDof_dt2 ... 2nd time derivative
    void ExtractDofValues(StructureOutputBlockVector& rDof_dt0, StructureOutputBlockVector& rDof_dt1,
                          StructureOutputBlockVector& rDof_dt2) const;

    //! @brief calculates the norm of the residual, can include weighting
    //! @param rResidual ... residual
    double CalculateNorm(const BlockFullVector<double>& rResidual) const;

    //! @brief postprocess (nodal dofs etc. and visualize a vtk file)
    //! @param rOutOfBalance ... out of balance values of the independent dofs (for disp dofs, this is the out of
    //! balance force)
    //! @remark rOutOfBalance here means Residual = ExternalForces - InternalForces
    void PostProcess(const StructureOutputBlockVector& rOutOfBalance);

    //! @brief returns the name of the restart file
    std::string GetRestartFileName() const;

    //! @brief sets the  time step for the time integration procedure (initial value)
    void SetTimeStep(double rTimeStep)
    {
        mTimeStep = rTimeStep;
    }

    //! @brief returns the  time step for the time integration procedure (current value)
    double GetTimeStep() const
    {
        return mTimeStep;
    }
    //! @brief sets the maximum time step for the time integration procedure
    void SetMaxTimeStep(double rMaxTimeStep)
    {
        mMaxTimeStep = rMaxTimeStep;
    }

    //! @brief returns the accumulated global time
    double GetTime() const
    {
        return mTime;
    }
    //! @brief sets the accumulated global time
    void SetTime(double rTime)
    {
        mTime = rTime;
    }

    //! @brief returns the maximum time step for the time integration procedure
    double GetMaxTimeStep() const
    {
        return mMaxTimeStep;
    }

    //! @brief sets the minimum time step for the time integration procedure
    void SetMinTimeStep(double rMinTimeStep)
    {
        mMinTimeStep = rMinTimeStep;
    }

    //! @brief sets the minimum time step for the time integration procedure
    void SetMinTimeStepPlot(double rMinTimeStepPlot)
    {
        mMinTimeStepPlot = rMinTimeStepPlot;
    }

    //! @brief returns the g
    std::vector<int> GetVecGroupNodesReactionForces() const
    {
        return mVecGroupNodesReactionForces;
    }

    //! @brief monitor the accelerations of a node
    //! @param rNodeId id of the node
    //! @param rResultId string identifying the result, this is used for the output file
    //! @return id of the result, so that it could be modified afterwards
    int AddResultNodeAccelerations(const std::string& rResultStr, int rNodeId);

    //! @brief monitor the displacements of a node
    //! @param rNodeId id of the node
    //! @param rResultId string identifying the result, this is used for the output file
    //! @return id of the result, so that it could be modified afterwards
    int AddResultNodeDisplacements(const std::string& rResultStr, int rNodeId);

    //! @brief monitor the time
    //! @param rResultId string identifying the result, this is used for the output file
    //! @param rGroupNodeId group id of the node group, for which the reaction forces (out of balance forces) should be
    //! calculated
    //! @return id of the result, so that it could be modified afterwards
    int AddResultGroupNodeForce(const std::string& rResultStr, int rGroupNodeId);

    //! @brief monitor the time
    //! @param rResultId string identifying the result, this is used for the output file
    //! @return id of the result, so that it could be modified afterwards
    int AddResultTime(const std::string& rResultStr);

    //! @brief monitor the integration point values in an element
    //! @param rResultId string identifying the result, this is used for the output file
    //! @param rElementId id of the element to be monitored
    //! @return id of the result, so that it could be modified afterwards
    int AddResultElementIpData(const std::string& rResultStr, int rElementId,
                               NuTo::IpData::eIpStaticDataType rIpDataType);

    //! @brief sets the result directory
    //! @param if delete is set, all the content of the directory will be removed
    void SetResultDirectory(std::string rResultDir, bool rDelete = false);

    //! @brief getter for the result directory
    std::string GetResultDirectory() const
    {
        return mResultDir;
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

protected:
    void ExportVisualizationFiles(const std::string& rResultDir, double rTime, int timeStep);

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

    double mTime = 0.; //!< accumulated time (in case several loadings are looked at, one after another)
    bool mAutomaticTimeStepping = false; //!< adapt the time step based on the number of iterations required (or
                                         //!decrease, if no convergence can be achieved)
    double mTimeStep = 0; //!< initial time step (or at the end this is the last time step used)

    double mMaxTimeStep = 1.;
    double mMinTimeStep = 0.;

    bool mMergeActiveDofValuesOrder1 = true; //!< if set to true, store velocities at the nodes in each time step
                                             //!(required when postprocessing velocities)
    bool mMergeActiveDofValuesOrder2 = false; //!< if set to true, store acceleration at the nodes in each time step
                                              //!(required when postprocessing accelerations)
    bool mCheckCoefficientMatrix = false; //!< if set to true, checks the coefficient matrix in each sub step
    bool mExportDataFileNodes = true; //!< if set to true, exports a data file for the nodes

    BlockScalar mToleranceResidual;
    //************************
    //* PostProcessing Stuff *
    //************************
    std::string mResultDir; //!< result directory

    boost::ptr_map<int, ResultBase> mResultMap; //!< specifies what to plot (displacements, reaction forces, etc.)

    int mLoadStep = 0; //!< load step number is increased after each converged step (used for successive output)
    int mTimeStepResult = 0; //!< time step number is increased each time a value is added to the result matrices
    int mTimeStepVTK = 0; //!< time step number is increased each time a vtk file is extracted
    double mMinTimeStepPlot = 0; //!< if the time between the current time step and the previous plotted step is larger
                                 //!than mMaxDeltaTimeStepPlot a vtk plot is performed
    double mLastTimePlot = -1e99; //!< last time when a vtk file was plotted
    int mIterationCount = 0; //!< iteration count


    std::vector<int> mVecGroupNodesReactionForces; //!< vector of groups of nodes for which the residual (corresponding
                                                   //!to the reaction forces induced by constraints) is given as output

    CallbackInterface* mCallback;

    std::vector<std::set<Node::eDof>> mStepActiveDofs; //!< stores wich Dofs are active in which calculation step

    bool mShowTime = true;
};
} // namespace NuTo
