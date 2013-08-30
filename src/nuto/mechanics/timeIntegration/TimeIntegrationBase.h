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
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"

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
    TimeIntegrationBase();

    //! @brief deconstructor
    virtual ~TimeIntegrationBase()
    {}

    //! @brief perform the time integration
    //! @param rStructure ... structure
    //! @param rTimeDelta ... length of the simulation
    virtual NuTo::Error::eError Solve(StructureBase& rStructure, double rTimeDelta)=0;

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
    double CalculateTimeDependentConstraintFactor(double rTimeDelta);

    //! @brief calculate the external force as a function of time delta
    //! @param curTime ... current time in the load step
    //! @param rLoad_j ... external load vector for the independent dofs
    //! @param rLoad_k ... external load vector for the dependent dofs
    void CalculateExternalLoad(StructureBase& rStructure, double curTime, NuTo::FullVector<double,Eigen::Dynamic>& rLoad_j, NuTo::FullVector<double,Eigen::Dynamic>& rLoad_k);

    //! @brief sets the nodes, for which displacements are to be monitored
    void CalculateOutputDispNodesPtr(StructureBase& rStructure);

    //! @brief calculate the external force as a function of time delta
    //! @ param rStructure ... structure
    //! @ param rPlotVector... data to be plotted, is append to the matrix and written to a file
    void PostProcess(StructureBase& rStructure, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rPlotVector);

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
    void SetGroupNodesReactionForces(NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> rVecGroupNodesReactionForces);

    //! @brief sets the nodes, for which displacements are to be monitored
    void SetOutputDispNodes(NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> rVecOutputDispNodes);

    //! @brief sets the minimum time step for the time integration procedure
    void SetPlotElementGroups(NuTo::FullVector<int,Eigen::Dynamic> rPlotElementGroups);

    //! @brief returns the g
    const NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> GetVecGroupNodesReactionForces()const
    {
    	return mVecGroupNodesReactionForces;
    }

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
    bool SetAutomaticTimeStepping()const
    {
    	return mAutomaticTimeStepping;
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
    //maximum time step (that's what we start with)
	double mMaxTimeStep;
    //minimum time step (for adaptive simulations)
	double mMinTimeStep;
	//load step number is increased after each converged step (used for successive output)
    int mLoadStep;
	//if the time between the current time step and the previous plotted step is larger than mMaxDeltaTimeStepPlot a vtk plot is performed
    double mMinTimeStepPlot;
    //last time when a vtk file was plotted
    double mLastTimePlot;
    //output data for all load steps
    FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mPlotMatrixAllLoadSteps;
    //output data for load steps that have a minimum distance of mMinTimeStepPlot
    FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> mPlotMatrixSelectedLoadSteps;
    //output data for load steps that have a minimum distance of mMinTimeStepPlot
    FullVector<int,Eigen::Dynamic> mPlotElementGroups;
    //result directory
    std::string mResultDir;
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
