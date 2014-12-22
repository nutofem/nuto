// $Id: TimeIntegrationBase.h 593 2012-02-03 23:50:23Z unger3 $

#ifndef TIMEINTEGRATIONBASE_H
#define TIMEINTEGRATIONBASE_H

#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <boost/ptr_container/ptr_map.hpp>

#include "nuto/base/NuToObject.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/mechanics/timeIntegration/ResultBase.h"

namespace NuTo
{
class StructureBase;
class NodeBase;
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
    virtual ~TimeIntegrationBase()
    {}

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    virtual NuTo::Error::eError Solve(double rTimeDelta)=0;

    //! @brief sets the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
    void ResetForNextLoad();

    //! @brief sets the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
    //! @param rTimeDependentConstraint ... constraint, whose rhs is increased as a function of time
    //! @param mTimeDependentConstraintFactor ... first row time, rhs of the constraint (linear interpolation in between afterwards linear extrapolation)
    void SetTimeDependentConstraint(int rTimeDependentConstraint, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& mTimeDependentConstraintFactor);

    //! @brief sets a scalar time dependent multiplication factor for the external loads
    //! @param rLoadRHSFactor ... first row time, second row scalar factor to calculate the external load (linear interpolation in between,  afterwards linear extrapolation)
    void SetTimeDependentLoadCase(int rLoadCase, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rLoadRHSFactor);

    //! @brief apply calculate the new rhs of the constraints as a function of the current time delta
    virtual double CalculateTimeDependentConstraintFactor(double rTimeDelta);

    //! @brief calculate the external force as a function of time delta
    //! @param curTime ... current time in the load step
    //! @param rLoad_j ... external load vector for the independent dofs
    //! @param rLoad_k ... external load vector for the dependent dofs
    void CalculateExternalLoad(StructureBase& rStructure, double curTime, NuTo::FullVector<double,Eigen::Dynamic>& rLoad_j, NuTo::FullVector<double,Eigen::Dynamic>& rLoad_k);

    //! @brief sets the nodes, for which displacements are to be monitored
    void CalculateOutputDispNodesPtr(StructureBase& rStructure);

    //! @brief postprocess (nodal dofs etc. and visualize a vtk file)
    //! @param rOutOfBalance_j ... out of balance values of the independent dofs (for disp dofs, this is the out of balance force)
    //! @param rOutOfBalance_k ... residual of the  dependent dofs
    void PostProcess(const FullVector<double,Eigen::Dynamic>& rOutOfBalance_j,
    		         const FullVector<double,Eigen::Dynamic>& rOutOfBalance_k);

    //! @brief visualizes the residual by merging it to the nodes, exporting the vtu file and restoring the inital state
    //! @param rResidual_j ... residual vector
    //! @param rActiveDofValues ... vector of active dof values
    //! @param rTimeStep ... time step number for file names
    void VisualizeResidual(
            const FullVector<double,Eigen::Dynamic>& rResidual_j,
            const FullVector<double, Eigen::Dynamic>& rActiveDofValues, int rTimeStep);

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

    //! @brief sets the minimum time step for the time integration procedure
    void SetPlotElementGroups(NuTo::FullVector<int,Eigen::Dynamic> rPlotElementGroups);

    //! @brief returns the g
    const NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> GetVecGroupNodesReactionForces()const
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

    //! @brief monitor the stress of all ips in an element (always 6 components times number of integration points)
    //! @param rResultId string identifying the result, this is used for the output file
    //! @param rElementId id of the element to be monitored
    //! @return id of the result, so that it could be modified afterwards
    int AddResultElementIpStress(const std::string& rResultStr, int rElementId);

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

    //! @brief sets the visualize residual flag (on or off)
    void SetVisualizeResidual(bool rVisualizeResidual)
    {
        mVisualizeResidual = rVisualizeResidual;
    }

    //! @brief returns if the visualize residual flag is turned on
    bool GetVisualizeResidual() const
    {
        return mVisualizeResidual;
    }

    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    virtual bool HasCriticalTimeStep()const = 0;

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    virtual double CalculateCriticalTimeStep()const = 0;

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
    TimeIntegrationBase(){};

    void ExportVisualizationFiles(const std::string& rResultDir, double rTime, int timeStep);


    //structure belonging to the time integration scheme
    StructureBase* mStructure;
    //constraint for displacement control (as a function of time)
    int mTimeDependentConstraint;
    //includes for each time step the rhs of the constraint mConstraintLoad
    //the time step is given by mTimeDelta/(mConstraintRHS.Rows()-1)
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mTimeDependentConstraintFactor;
    //time dependent load case number
    int mTimeDependentLoadCase;
	//includes for each time step the scalar factor for the load case
    //the time step is given relative to mTimeDelta
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mTimeDependentLoadFactor;
    //external load vectors (static and time dependent)
	NuTo::FullVector<double,Eigen::Dynamic> mLoadVectorStatic_j,mLoadVectorStatic_k,mLoadVectorTimeDependent_j,mLoadVectorTimeDependent_k;
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

    //if set to true, prints an additional result files with residual values
    bool mVisualizeResidual;

    //groups of elements to be plotted separately
    FullVector<int,Eigen::Dynamic> mPlotElementGroups;

    // vector of groups of nodes for which the residual (corresponding to the reaction forces induced by constraints) is given as output
    NuTo::FullVector<int,Eigen::Dynamic> mVecGroupNodesReactionForces;
    // vector of nodes for the output of the dofs
    std::vector<NodeBase*> mVecOutputDispNodesPtr;
    NuTo::FullVector<int,Eigen::Dynamic> mVecOutputDispNodesInt;


};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::TimeIntegrationBase)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // TIMEINTEGRATIONBASE_H
