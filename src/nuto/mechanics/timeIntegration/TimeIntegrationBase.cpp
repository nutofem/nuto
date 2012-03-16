// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#include <boost/ptr_container/serialize_ptr_list.hpp>
#endif // ENABLE_SERIALIZATION

#include "boost/filesystem.hpp"
#include <iostream>
#include <fstream>
#ifdef SHOW_TIME
#include <ctime>
#endif

# ifdef _OPENMP
#include <omp.h>
# endif


#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/StructureBase.h"


NuTo::TimeIntegrationBase::TimeIntegrationBase()  : NuTo::NuToObject::NuToObject()
{
	mTime = 0.;
	mLoadStep = 1;
	mMaxTimeStep = 1;
	mMinTimeStep = 0;
    mLoadStep = 0;
    mMinTimeStepPlot = 0;
    mLastTimePlot = 0;
    mAutomaticTimeStepping = false;
 	ResetForNextLoad();
}

//! @brief sets the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
void NuTo::TimeIntegrationBase::ResetForNextLoad()
{
	mConstraintLoad = -1;
	mConstraintRHS.Resize(0,0);
	mLoadRHSFactor.Resize(0,0);
}

//! @brief sets the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
//! @param rConstraintLoad ... constraint, whose rhs is increased as a function of time
//! @param rConstraintRHS ... first row time, rhs of the constraint (linear interpolation in between afterwards, constant)
void NuTo::TimeIntegrationBase::SetDisplacements(int rConstraintLoad, const NuTo::FullMatrix<double>& rConstraintRHS)
{
	if (rConstraintRHS.GetNumColumns()!=2)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetDisplacements] number of columns must be 2, first column contains the time, second column contains the corresponding value.");
	if (rConstraintRHS.GetNumRows()<2)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetDisplacements] number of rows must be at least 2.");
	if (rConstraintRHS(0,0)!=0)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetDisplacements] the first time should always be zero.");

	//check, if the time is monotonically increasing
	for (int count=0; count<rConstraintRHS.GetNumRows()-1; count++)
	{
		if (rConstraintRHS(count,0)>=rConstraintRHS(count+1,0))
			throw MechanicsException("[NuTo::TimeIntegrationBase::SetDisplacements] time has to increase monotonically.");
	}

	mConstraintLoad = rConstraintLoad;
	mConstraintRHS = rConstraintRHS;
}

//! @brief sets a scalar time dependent multiplication factor for the external loads
//! @param rLoadRHSFactor ... first row time, second row scalar factor to calculate the external load (linear interpolation in between, afterwards constant)
void NuTo::TimeIntegrationBase::SetExternalLoads(const NuTo::FullMatrix<double>& rLoadRHSFactor)
{
	if (rLoadRHSFactor.GetNumColumns()!=2)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetExternalLoads] number of columns must be 2, first column contains the time, second column contains the corresponding value.");
	if (rLoadRHSFactor.GetNumRows()<2)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetExternalLoads] number of rows must be at least 2.");
	if (rLoadRHSFactor(0,0)!=0)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetExternalLoads] the first time should always be zero.");

	//check, if the time is monotonically increasing
	for (int count=0; count<rLoadRHSFactor.GetNumRows()-1; count++)
	{
		if (rLoadRHSFactor(count,0)>=rLoadRHSFactor(count+1,0))
			throw MechanicsException("[NuTo::TimeIntegrationBase::SetExternalLoads] time has to increase monotonically.");
	}

	mLoadRHSFactor = rLoadRHSFactor;
}

//! @brief apply the new rhs of the constraints as a function of the current time delta
void NuTo::TimeIntegrationBase::ConstraintsCalculateRHSAndApply(StructureBase& rStructure, double curTime)
{
    //calculate the two corresponding time steps between which a linear interpolation is performed
	if (mConstraintRHS.GetNumRows()!=0)
	{
        int curStep(0);
        while (mConstraintRHS(curStep,0)<curTime && curStep<mConstraintRHS.GetNumRows()-1)
        	curStep++;
		if (curStep==(mConstraintRHS.GetNumRows()-1))
			curStep--;

		//extract the two data points
		double s1 = mConstraintRHS(curStep,1);
		double s2 = mConstraintRHS(curStep+1,1);
		double t1 = mConstraintRHS(curStep,0);
		double t2 = mConstraintRHS(curStep+1,0);

		double s = 	s1 + (s2-s1)/(t2-t1) * (curTime-t1);

		rStructure.ConstraintSetRHS(mConstraintLoad,s);
	}
}

//! @brief calculate the external force as a function of time delta
//! @param curTime ... current time (within the loadstep)
//! @param rLoad_j ... external load vector for the independent dofs
//! @param rLoad_k ... external load vector for the dependent dofs
void NuTo::TimeIntegrationBase::CalculateExternalLoad(StructureBase& rStructure, double curTime, NuTo::FullMatrix<double>& rLoad_j, NuTo::FullMatrix<double>& rLoad_k)
{
	if (mLoadRHSFactor.GetNumRows()!=0)
	{
        int curStep(0);
        while (mLoadRHSFactor(curStep,0)<curTime && curStep<mLoadRHSFactor.GetNumRows()-1)
        	curStep++;
		if (curStep==(mLoadRHSFactor.GetNumRows()-1))
			curStep--;

		//extract the two data points
		double s1 = mLoadRHSFactor(curStep,1);
		double s2 = mLoadRHSFactor(curStep+1,1);
		double t1 = mLoadRHSFactor(curStep,0);
		double t2 = mLoadRHSFactor(curStep+1,0);

		double s = 	s1 + (s2-s1)/(t2-t1) * (curTime-t1);

		//in theory, this vector could be stored, but I guess the recalculation is in general not very expensive
		rStructure.BuildGlobalExternalLoadVector(rLoad_j,rLoad_k);
		rLoad_j*=s;
		rLoad_k*=s;
	}
	else
	{
		rStructure.BuildGlobalExternalLoadVector(rLoad_j,rLoad_k);
	}
}



//! @brief calculate the external force as a function of time delta
//! @ param rStructure ... structure
//! @ param rPlotVector... data to be plotted, is append to the matrix and written to a file
void NuTo::TimeIntegrationBase::PostProcess(StructureBase& rStructure, FullMatrix<double>& rPlotVector)
{
	if (mResultDir.length()>0)
	{
		if (mPlotMatrixAllLoadSteps.GetNumRows()==0)
			mPlotMatrixAllLoadSteps = rPlotVector;
		else
			mPlotMatrixAllLoadSteps.AppendRows(rPlotVector);
		boost::filesystem::path resultFile(mResultDir);
		resultFile /= std::string("resultAllLoadSteps.dat");
		mPlotMatrixAllLoadSteps.WriteToFile(resultFile.string(), std::string("  "));

		if (mTime-mLastTimePlot>=mMinTimeStepPlot)
		{
			if (mPlotMatrixSelectedLoadSteps.GetNumRows()==0)
				mPlotMatrixSelectedLoadSteps = rPlotVector;
			else
			    mPlotMatrixSelectedLoadSteps.AppendRows(rPlotVector);
			boost::filesystem::path resultFileMod(mResultDir);
			resultFileMod /= std::string("resultSelectedLoadSteps.dat");
			mPlotMatrixSelectedLoadSteps.WriteToFile(resultFileMod.string(), std::string("  "));

			int curModLoadStep=mPlotMatrixSelectedLoadSteps.GetNumRows();

#ifdef ENABLE_VISUALIZE
			//plot the solution vtk file
			std::stringstream ssLoadStep;
			ssLoadStep << curModLoadStep;
			resultFile = mResultDir;
			resultFile /= std::string("Elements")+ssLoadStep.str()+std::string(".vtu");
			rStructure.ExportVtkDataFileElements(resultFile.string(),true);
			resultFile = mResultDir;
			resultFile /= std::string("Nodes")+ssLoadStep.str()+std::string(".vtu");
			rStructure.ExportVtkDataFileNodes(resultFile.string(),true);

			//write an additional pvd file
			resultFile = mResultDir;
			resultFile /= std::string("Elements.pvd");
		    std::fstream file;
		    if (curModLoadStep==1)
		    {
		    	file.open(resultFile.c_str(), std::fstream::out);
		    }
		    else
		    {
		    	file.open(resultFile.c_str(), std::fstream::out | std::fstream::in |std::ios_base::ate);
		    }

		    if (!file.is_open())
		    {
		    	throw NuTo::MechanicsException(std::string("[NuTo::TimeIntegrationBase::PostProcess] Error opening file ")+resultFile.c_str());
		    }
		    std::stringstream endOfXML;
		    endOfXML << "</Collection>" << std::endl;
		    endOfXML << "</VTKFile>" << std::endl;
		    if (curModLoadStep==1)
		    {
				// header /////////////////////////////////////////////////////////////////
				file << "<?xml version=\"1.0\"?>" << std::endl;
				file << "<VTKFile type=\"Collection\">" << std::endl;
				file << "<Collection>" << std::endl;
				file << "<DataSet timestep=\""<< mTime << "\" file=\"Elements" << curModLoadStep << ".vtu\"/>" << std::endl;
		    }
		    else
		    {
			    //delete the last part of the xml file
		    	file.seekp (-endOfXML.str().length(),std::ios_base::end);
				file << "<DataSet timestep=\""<< mTime << "\" file=\"Elements" << curModLoadStep << ".vtu\"/>" << std::endl;
		    }
		    file << endOfXML.str();
		    file.close();
#endif
			mLastTimePlot = mTime;
		}
    }
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::TimeIntegrationBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::TimeIntegrationBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::TimeIntegrationBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::TimeIntegrationBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::TimeIntegrationBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::TimeIntegrationBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::TimeIntegrationBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    mLogger << "start serialization of TimeIntegrationBase" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
       & BOOST_SERIALIZATION_NVP(mConstraintLoad)
       & BOOST_SERIALIZATION_NVP(mConstraintRHS)
       & BOOST_SERIALIZATION_NVP(mLoadRHSFactor)
       & BOOST_SERIALIZATION_NVP(mTime)
       & BOOST_SERIALIZATION_NVP(mLoadStep)
       & BOOST_SERIALIZATION_NVP(mMaxTimeStep)
       & BOOST_SERIALIZATION_NVP(mMinTimeStep)
       & BOOST_SERIALIZATION_NVP(mLoadStep)
       & BOOST_SERIALIZATION_NVP(mMinTimeStepPlot)
       & BOOST_SERIALIZATION_NVP(mLastTimePlot)
       & BOOST_SERIALIZATION_NVP(mPlotMatrixAllLoadSteps)
       & BOOST_SERIALIZATION_NVP(mPlotMatrixSelectedLoadSteps)
       & BOOST_SERIALIZATION_NVP(mResultDir)
       & BOOST_SERIALIZATION_NVP(mAutomaticTimeStepping);
#ifdef DEBUG_SERIALIZATION
    mLogger << "finish serialization of structure base" << "\n";
#endif
}
#endif  // ENABLE_SERIALIZATION

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::TimeIntegrationBase::Info()const
{
}

//! @brief sets the result directory
void NuTo::TimeIntegrationBase::SetResultDirectory(std::string rResultDir)
{
	mResultDir = rResultDir;
    //delete result directory
    if (boost::filesystem::exists(rResultDir))    // does p actually exist?
    {
        if (boost::filesystem::is_directory(rResultDir))      // is p a directory?
        {
        	boost::filesystem::remove_all(rResultDir);
        }
    }
    // create result directory
    boost::filesystem::create_directory(rResultDir);
}

//! @brief sets the minimum time step for the time integration procedure
void NuTo::TimeIntegrationBase::SetGroupNodesReactionForces(NuTo::FullMatrix<int> rVecGroupNodesReactionForces)
{
	if (rVecGroupNodesReactionForces.GetNumColumns()!=1)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetGroupNodesReactionForces] vector (must have a single column).");
	if (rVecGroupNodesReactionForces.GetNumRows()<1)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetGroupNodesReactionForces] vector must have at least a single row.");
	mVecGroupNodesReactionForces = rVecGroupNodesReactionForces;
}

