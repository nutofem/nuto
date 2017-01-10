// $Id: TimeIntegrationBase.h 593 2012-02-03 23:50:23Z unger3 $

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <vector>
#include <boost/ptr_container/ptr_map.hpp>


// parent class
#include "base/NuToObject.h"

// member
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

namespace NuTo
{
class CallbackInterface;
class NodeBase;
class ResultBase;
class StructureBase;
class TimeDependencyBase;
enum class eError;

template <typename T> class BlockFullVector;

namespace IpData
{
    enum class eIpStaticDataType;
}// namespace IpData

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
    TimeIntegrationBase(StructureBase* rStructure);

    //! @brief deconstructor
    virtual ~TimeIntegrationBase();

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    virtual NuTo::eError Solve(double rTimeDelta)=0;

    //! @brief sets the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
    void ResetForNextLoad();

    //! @brief Adds the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
    //! @param rTimeDependentConstraint ... constraint, whose rhs is increased as a function of time
    //! @param mTimeDependentConstraintFactor ... first row time, rhs of the constraint (linear interpolation in between afterwards linear extrapolation)
    void AddTimeDependentConstraint(int rTimeDependentConstraint, const Eigen::MatrixXd& mTimeDependentConstraintFactor);

    //! @brief Adds the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
    //! @param rTimeDependentConstraint ... constraint, whose rhs is increased as a function of time
    //! @param rTimeDependentConstraintFunction ... function that calculates the time dependent constraint factor for the current time step
    void AddTimeDependentConstraintFunction(int rTimeDependentConstraint, const std::function<double (double rTime)>& rTimeDependentConstraintFunction);


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
    //! @param rLoadRHSFactor ... first row time, second row scalar factor to calculate the external load (linear interpolation in between,  afterwards linear extrapolation)
    void SetTimeDependentLoadCase(int rLoadCase, const Eigen::MatrixXd& rLoadRHSFactor);

    //! @brief apply calculate the new rhs of the constraints as a function of the current time delta
    virtual double CalculateTimeDependentConstraintFactor(double rTimeDelta);

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
    void ExtractDofValues(StructureOutputBlockVector& rDof_dt0, StructureOutputBlockVector& rDof_dt1, StructureOutputBlockVector& rDof_dt2) const;

    //! @brief calculates the norm of the residual, can include weighting
    //! @param rResidual ... residual
    double CalculateNorm(const BlockFullVector<double>& rResidual) const;

    //! @brief postprocess (nodal dofs etc. and visualize a vtk file)
    //! @param rOutOfBalance ... out of balance values of the independent dofs (for disp dofs, this is the out of balance force)
    void PostProcess(const StructureOutputBlockVector& rOutOfBalance);

    //! @brief sets the  time step for the time integration procedure (initial value)
    void SetTimeStep(double rTimeStep)
    {
        mTimeStep = rTimeStep;
    }

    //! @brief returns the  time step for the time integration procedure (current value)
    double GetTimeStep()const
    {
        return mTimeStep;
    }

    //! @brief sets the maximum time step for the time integration procedure
    void SetMaxTimeStep(double rMaxTimeStep)
    {
    	mMaxTimeStep = rMaxTimeStep;
    }

    //! @brief returns the maximum time step for the time integration procedure
    double GetMaxTimeStep()const
    {
    	return mMaxTimeStep;
    }

    //! @brief sets the minimum time step for the time integration procedure
    void SetMinTimeStep(double rMinTimeStep)
    {
    	mMinTimeStep = rMinTimeStep;
    }

    //! @brief returns the minimum time step for the time integration procedure
    double GetMinTimeStep()const
    {
    	return mMinTimeStep;
    }

    //! @brief sets the minimum time step for the time integration procedure
    void SetMinTimeStepPlot(double rMinTimeStepPlot)
    {
    	mMinTimeStepPlot = rMinTimeStepPlot;
    }

    //! @brief returns the minimum time step for the time integration procedure
    double GetMinTimeStepPlot()const
    {
    	return mMinTimeStepPlot;
    }

    void SetPlotElementGroups(std::vector<int> rPlotElementGroups);

    //! @brief returns the g
    std::vector<int> GetVecGroupNodesReactionForces()const
    {
    	return mVecGroupNodesReactionForces;
    }

    //! @brief monitor the accelerations of a node
    //! @param rNodeId id of the node
    //! @param rResultId string identifying the result, this is used for the output file
    //! @return id of the result, so that it could be modified afterwards
    int AddResultNodeAccelerations(const std::string& rResultStr, int rNodeId );

    //! @brief monitor the displacements of a node
    //! @param rNodeId id of the node
    //! @param rResultId string identifying the result, this is used for the output file
    //! @return id of the result, so that it could be modified afterwards
    int AddResultNodeDisplacements(const std::string& rResultStr,int rNodeId);

    //! @brief monitor the time
    //! @param rResultId string identifying the result, this is used for the output file
    //! @param rGroupNodeId group id of the node group, for which the reaction forces (out of balance forces) should be calculated
    //! @return id of the result, so that it could be modified afterwards
    int AddResultGroupNodeForce(const std::string& rResultStr,int rGroupNodeId);

    //! @brief monitor the time
    //! @param rResultId string identifying the result, this is used for the output file
    //! @return id of the result, so that it could be modified afterwards
    int AddResultTime(const std::string& rResultStr);

    //! @brief monitor the integration point values in an element
    //! @param rResultId string identifying the result, this is used for the output file
    //! @param rElementId id of the element to be monitored
    //! @return id of the result, so that it could be modified afterwards
    int AddResultElementIpData(const std::string& rResultStr, int rElementId, NuTo::IpData::eIpStaticDataType rIpDataType);

    //! @brief sets the result directory
    //! @param if delete is set, all the content of the directory will be removed
    void SetResultDirectory(std::string rResultDir, bool rDelete);

    //! @brief returns the result directory
    std::string GetResultDirectory()const
    {
    	return mResultDir;
    }

    //! @brief sets automatic time stepping (on or off)
    void SetAutomaticTimeStepping(bool rAutomaticTimeStepping)
    {
    	mAutomaticTimeStepping = rAutomaticTimeStepping;
    }

    //! @brief returns if automatic time stepping is turned on
    bool GetAutomaticTimeStepping()const
    {
    	return mAutomaticTimeStepping;
    }

    //! @brief sets the coefficient matrix check (on or off)
    void SetCheckCoefficientMatrix(bool rCheckCoefficientMatrix)
    {
        mCheckCoefficientMatrix = rCheckCoefficientMatrix;
    }

    //! @brief returns if the coefficient matrix check is turned on
    bool GetCheckCoefficientMatrix() const
    {
        return mCheckCoefficientMatrix;
    }

    //! @brief Sets the export data file nodes bool
    void SetExportDataFileNodes(bool rExportDataFileNodes)
    {
        mExportDataFileNodes = rExportDataFileNodes;
    }

    //! @brief Gets the export data file nodes bool
    bool GetExportDataFileNodes() const
    {
        return mExportDataFileNodes;
    }

    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    virtual bool HasCriticalTimeStep()const = 0;

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    virtual double CalculateCriticalTimeStep()const = 0;

    void ConnectCallback(CallbackInterface* rCallback)
    {
        mCallback = rCallback;
    }

    //! @brief ... Adds a calculation step to each timestep
    //! param rActiveDofs ... active Dofs of the calculation step
    void AddCalculationStep(const std::set<Node::eDof> &rActiveDofs);

    //! @brief ... Sets the number of calculation steps per timestep
    //! param rNumSteps ... number of calculation steps per timestep
    void SetNumCalculationSteps(int rNumSteps);

    //! @brief ... Sets the active Dofs of a calculation step
    //! param rStepNum ... step number
    //! param rActiveDofs ... active Dofs of the calculation step
    void SetActiveDofsCalculationStep(int rStepNum, const std::set<Node::eDof> &rActiveDofs);

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
    //empty private construct required for serialization

#ifdef ENABLE_SERIALIZATION
    TimeIntegrationBase() : mLoadVectorStatic(DofStatus()), mLoadVectorTimeDependent(DofStatus()) {};

    //empty private construct required for serialization
    NewmarkDirect() : mToleranceResidual(DofStatus()) {};
#endif // ENABLE_SERIALIZATION
    void ExportVisualizationFiles(const std::string& rResultDir, double rTime, int timeStep);

    const BlockFullVector<double>& UpdateAndGetConstraintRHS(double rCurrentTime);

    const BlockFullVector<double>& UpdateAndGetAndMergeConstraintRHS(double rCurrentTime, StructureOutputBlockVector& rDof_dt0);

    //structure belonging to the time integration scheme
    StructureBase* mStructure;


    // DEPRECATED BLOCK BEGIN

        //constraint for displacement control (as a function of time)
        int mTimeDependentConstraint;
        //includes for each time step the rhs of the constraint mConstraintLoad
        //the time step is given by mTimeDelta/(mConstraintRHS.Rows()-1)
        Eigen::MatrixXd mTimeDependentConstraintFactor;

    // DEPRECATED BLOCK END

    // Saves the IDs of all time dependent constraints
    std::map<int,std::shared_ptr<TimeDependencyBase>> mMapTimeDependentConstraint;



    //time dependent load case number
    int mTimeDependentLoadCase;
	//includes for each time step the scalar factor for the load case
    //the time step is given relative to mTimeDelta
	Eigen::MatrixXd mTimeDependentLoadFactor;
    //external load vectors (static and time dependent)
	NuTo::StructureOutputBlockVector mLoadVectorStatic;
	NuTo::StructureOutputBlockVector mLoadVectorTimeDependent;

	//accumulated time (in case several loadings are looked at, one after another)
	double mTime;
    //adapt the time step based on the number of iterations required (or decrease, if no convergence can be achieved)
	bool mAutomaticTimeStepping;
	//initial time step (or at the end this is the last time step used)
	double mTimeStep;
    //maximum time step (for adaptive simulations)
	double mMaxTimeStep;
    //minimum time step (for adaptive simulations)
	double mMinTimeStep;
	//if set to true, store velocities at the nodes in each time step (required when postprocessing velocities)
	bool mMergeActiveDofValuesOrder1;
	//if set to true, store acceleration at the nodes in each time step (required when postprocessing accelerations)
	bool mMergeActiveDofValuesOrder2;
	//if set to true, checks the coefficient matrix in each sub step
	bool mCheckCoefficientMatrix;
    //! @brief If set to true, exports a data file for the nodes
    bool mExportDataFileNodes = true;

    BlockScalar mToleranceResidual;
	//************************
	//* PostProcessing Stuff *
	//************************
    //result directory
    std::string mResultDir;

	//specifies what to plot (displacements, reaction forces, etc.)
    boost::ptr_map<int,ResultBase> mResultMap;

	//load step number is increased after each converged step (used for successive output)
    int mLoadStep;
    //time step number is increased each time a value is added to the result matrices
    int mTimeStepResult;
    //time step number is increased each time a vtk file is extracted
    int mTimeStepVTK;
	//if the time between the current time step and the previous plotted step is larger than mMaxDeltaTimeStepPlot a vtk plot is performed
    double mMinTimeStepPlot;
    //last time when a vtk file was plotted
    double mLastTimePlot;

    //groups of elements to be plotted separately
    std::vector<int> mPlotElementGroups;

    // vector of groups of nodes for which the residual (corresponding to the reaction forces induced by constraints) is given as output
    std::vector<int> mVecGroupNodesReactionForces;

    CallbackInterface* mCallback;

    //! @brief Stores wich Dofs are active in which calculation step
    std::vector<std::set<Node::eDof>> mStepActiveDofs;

};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::TimeIntegrationBase)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
