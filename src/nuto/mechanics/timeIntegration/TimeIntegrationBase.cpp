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
 	mTimeDependentConstraint = -1;
 	mTimeDependentLoadCase = -1;
 	ResetForNextLoad();
}

//! @brief sets the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
void NuTo::TimeIntegrationBase::ResetForNextLoad()
{
	mTimeDependentConstraint = -1;
	mTimeDependentConstraintFactor.Resize(0,0);
 	mTimeDependentLoadCase = -1;
 	mTimeDependentLoadFactor.Resize(0,0);
}

//! @brief sets the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
//! @param rTimeDependentConstraint ... constraint, whose rhs is increased as a function of time
//! @param mTimeDependentConstraintFactor ... first row time, rhs of the constraint (linear interpolation in between afterwards, constant)
void NuTo::TimeIntegrationBase::SetTimeDependentConstraint(int rTimeDependentConstraint, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rTimeDependentConstraintFactor)
{
	if (rTimeDependentConstraintFactor.GetNumColumns()!=2)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetDisplacements] number of columns must be 2, first column contains the time, second column contains the corresponding rhs.");
	if (rTimeDependentConstraintFactor.GetNumRows()<2)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetDisplacements] number of rows must be at least 2.");
	if (rTimeDependentConstraintFactor(0,0)!=0)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetDisplacements] the first time should always be zero.");

	//check, if the time is monotonically increasing
	for (int count=0; count<rTimeDependentConstraintFactor.GetNumRows()-1; count++)
	{
		if (rTimeDependentConstraintFactor(count,0)>=rTimeDependentConstraintFactor(count+1,0))
			throw MechanicsException("[NuTo::TimeIntegrationBase::SetDisplacements] time has to increase monotonically.");
	}

	mTimeDependentConstraint = rTimeDependentConstraint;
	mTimeDependentConstraintFactor = rTimeDependentConstraintFactor;
}

//! @brief sets a scalar time dependent multiplication factor for the external loads
//! @param rTimeDependentLoadFactor ... first row time, second row scalar factor to calculate the external load (linear interpolation in between, afterwards constant)
void NuTo::TimeIntegrationBase::SetTimeDependentLoadCase(int rTimeDependentLoadCase, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rTimeDependentLoadFactor)
{
	if (rTimeDependentLoadFactor.GetNumColumns()!=2)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetExternalLoads] number of columns must be 2, first column contains the time, second column contains the corresponding value.");
	if (rTimeDependentLoadFactor.GetNumRows()<2)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetExternalLoads] number of rows must be at least 2.");
	if (rTimeDependentLoadFactor(0,0)!=0)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetExternalLoads] the first time should always be zero.");

	//check, if the time is monotonically increasing
	for (int count=0; count<rTimeDependentLoadFactor.GetNumRows()-1; count++)
	{
		if (rTimeDependentLoadFactor(count,0)>=rTimeDependentLoadFactor(count+1,0))
			throw MechanicsException("[NuTo::TimeIntegrationBase::SetExternalLoads] time has to increase monotonically.");
	}

	mTimeDependentLoadFactor = rTimeDependentLoadFactor;
	mTimeDependentLoadCase = rTimeDependentLoadCase;
}

//! @brief apply the new rhs of the constraints as a function of the current time delta
double NuTo::TimeIntegrationBase::CalculateTimeDependentConstraintFactor(double curTime)
{
    //calculate the two corresponding time steps between which a linear interpolation is performed
	if (mTimeDependentConstraintFactor.GetNumRows()!=0)
	{
        int curStep(0);
        while (mTimeDependentConstraintFactor(curStep,0)<curTime && curStep<mTimeDependentConstraintFactor.GetNumRows()-1)
        	curStep++;
		if (curStep==0)
			curStep++;

		//extract the two data points
		double s1 = mTimeDependentConstraintFactor(curStep-1,1);
		double s2 = mTimeDependentConstraintFactor(curStep,1);
		double t1 = mTimeDependentConstraintFactor(curStep-1,0);
		double t2 = mTimeDependentConstraintFactor(curStep,0);

		return s1 + (s2-s1)/(t2-t1) * (curTime-t1);
	}
	return 0;
}

//! @brief calculate the external force as a function of time delta
//! @param curTime ... current time (within the loadstep)
//! @param rLoad_j ... external load vector for the independent dofs
//! @param rLoad_k ... external load vector for the dependent dofs
void NuTo::TimeIntegrationBase::CalculateExternalLoad(StructureBase& rStructure, double curTime, NuTo::FullVector<double,Eigen::Dynamic>& rLoad_j, NuTo::FullVector<double,Eigen::Dynamic>& rLoad_k)
{
	if (mLoadVectorStatic_j.GetNumRows()==0 && mLoadVectorStatic_k.GetNumRows()==0 &&
		mLoadVectorTimeDependent_j.GetNumRows()==0 && mLoadVectorTimeDependent_k.GetNumRows()==0)
	{
		//build static and dynamic load vector
		NuTo::FullVector<double,Eigen::Dynamic> tmp_j,tmp_k;

		mLoadVectorStatic_j.Resize(rStructure.GetNumActiveDofs());
		mLoadVectorStatic_k.Resize(rStructure.GetNumDofs()-rStructure.GetNumActiveDofs());
		for (int count=0; count<rStructure.GetNumLoadCases(); count++)
		{
			rStructure.BuildGlobalExternalLoadVector(count,tmp_j,tmp_k);
			if (count==mTimeDependentLoadCase)
			{
				mLoadVectorTimeDependent_j=tmp_j;
				mLoadVectorTimeDependent_k=tmp_k;
			}
			else
			{
				mLoadVectorStatic_j+=tmp_j;
				mLoadVectorStatic_k+=tmp_k;
			}
			//std::cout << "tmp_j\n" << tmp_j <<  std::endl;
			//std::cout << "tmp_k\n" << tmp_k <<  std::endl;
			std::cout << "sum of loads for loadcase " << count << " is " << tmp_j.ColumnwiseSum() + tmp_k.ColumnwiseSum() << std::endl;
		}
	}

	if (mTimeDependentLoadCase!=-1)
	{
		if (mTimeDependentLoadFactor.GetNumRows()==0)
		{
			throw MechanicsException("[NuTo::TimeIntegrationBase::CalculateExternalLoad] TimeDependentLoadFactor not set.");
		}
        int curStep(0);
        while (mTimeDependentLoadFactor(curStep,0)<curTime && curStep<mTimeDependentLoadFactor.GetNumRows()-1)
        	curStep++;
		if (curStep==0)
			curStep++;

		//extract the two data points
		double s1 = mTimeDependentLoadFactor(curStep-1,1);
		double s2 = mTimeDependentLoadFactor(curStep,1);
		double t1 = mTimeDependentLoadFactor(curStep-1,0);
		double t2 = mTimeDependentLoadFactor(curStep,0);

		double s = 	s1 + (s2-s1)/(t2-t1) * (curTime-t1);

		rLoad_j=mLoadVectorStatic_j+mLoadVectorTimeDependent_j*s;
		rLoad_k=mLoadVectorStatic_k+mLoadVectorTimeDependent_k*s;
	}
	else
	{
		rLoad_j=mLoadVectorStatic_j;
		rLoad_k=mLoadVectorStatic_k;
	}
}



//! @brief calculate the external force as a function of time delta
//! @ param rStructure ... structure
//! @ param rPlotVector... data to be plotted, is append to the matrix and written to a file
void NuTo::TimeIntegrationBase::PostProcess(StructureBase& rStructure, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rPlotVector)
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
			resultFile /= std::string("Nodes")+ssLoadStep.str()+std::string(".vtu");
			rStructure.ExportVtkDataFileNodes(resultFile.string(),true);

			if (mPlotElementGroups.GetNumRows()==0)
			{
				//plot all elements
				resultFile = mResultDir;
				resultFile /= std::string("Elements")+ssLoadStep.str()+std::string(".vtu");
				rStructure.ExportVtkDataFileElements(resultFile.string(),true);

				//write an additional pvd file
				resultFile = mResultDir;
				resultFile /= std::string("Elements.pvd");
			    std::fstream file;
			    if (curModLoadStep==1)
			    {
			    	file.open(resultFile.string(), std::fstream::out);
			    }
			    else
			    {
			    	file.open(resultFile.string(), std::fstream::out | std::fstream::in |std::ios_base::ate);
			    }

			    if (!file.is_open())
			    {
			    	throw NuTo::MechanicsException(std::string("[NuTo::TimeIntegrationBase::PostProcess] Error opening file ")+resultFile.string());
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
			}
			else
			{
				//plot all groups separately
				for (int countGroupElement=0; countGroupElement<mPlotElementGroups.GetNumRows();countGroupElement++)
				{
					std::stringstream ssGroup;
					ssGroup << mPlotElementGroups(countGroupElement);

					//plot all elements
					resultFile = mResultDir;
					resultFile /= std::string("Group") + ssGroup.str() + std::string("_Elements")+ssLoadStep.str()+std::string(".vtu");
					rStructure.ElementGroupExportVtkDataFile(mPlotElementGroups(countGroupElement), resultFile.string(),true);

					//write an additional pvd file
					resultFile = mResultDir;
					resultFile /= std::string("Group") + ssGroup.str() + std::string("_ElementsAll")+std::string(".pvd");
				    std::fstream file;
				    if (curModLoadStep==1)
				    {
				    	file.open(resultFile.string(), std::fstream::out);
				    }
				    else
				    {
				    	file.open(resultFile.string(), std::fstream::out | std::fstream::in |std::ios_base::ate);
				    }

				    if (!file.is_open())
				    {
				    	throw NuTo::MechanicsException(std::string("[NuTo::TimeIntegrationBase::PostProcess] Error opening file ")+resultFile.string());
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
				    }
				    else
				    {
					    //delete the last part of the xml file
				    	file.seekp (-endOfXML.str().length(),std::ios_base::end);
				    }
					file << "<DataSet timestep=\""<< mTime << "\" file=\"Group" << mPlotElementGroups(countGroupElement) << "_Elements" << curModLoadStep << ".vtu\"/>" << std::endl;
				    file << endOfXML.str();
				    file.close();
				}
			}
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
    std::cout << "start serialization of TimeIntegrationBase" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
       & BOOST_SERIALIZATION_NVP(mTimeDependentConstraint)
       & BOOST_SERIALIZATION_NVP(mTimeDependentLoadFactor)
       & BOOST_SERIALIZATION_NVP(mTimeDependentLoadCase)
       & BOOST_SERIALIZATION_NVP(mTimeDependentLoadFactor)
       & BOOST_SERIALIZATION_NVP(mLoadVectorStatic_j)
       & BOOST_SERIALIZATION_NVP(mLoadVectorStatic_k)
       & BOOST_SERIALIZATION_NVP(mLoadVectorTimeDependent_j)
       & BOOST_SERIALIZATION_NVP(mLoadVectorTimeDependent_k)
       & BOOST_SERIALIZATION_NVP(mTime)
       & BOOST_SERIALIZATION_NVP(mLoadStep)
       & BOOST_SERIALIZATION_NVP(mMaxTimeStep)
       & BOOST_SERIALIZATION_NVP(mMinTimeStep)
       & BOOST_SERIALIZATION_NVP(mLoadStep)
       & BOOST_SERIALIZATION_NVP(mMinTimeStepPlot)
       & BOOST_SERIALIZATION_NVP(mLastTimePlot)
       & BOOST_SERIALIZATION_NVP(mPlotMatrixAllLoadSteps)
       & BOOST_SERIALIZATION_NVP(mPlotMatrixSelectedLoadSteps)
       & BOOST_SERIALIZATION_NVP(mPlotElementGroups)
       & BOOST_SERIALIZATION_NVP(mResultDir)
       & BOOST_SERIALIZATION_NVP(mAutomaticTimeStepping);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of structure base" << "\n";
#endif
}
#endif  // ENABLE_SERIALIZATION

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::TimeIntegrationBase::Info()const
{
}

//! @brief sets the result directory
//! @param if delete is set, all the content of the directory will be removed
void NuTo::TimeIntegrationBase::SetResultDirectory(std::string rResultDir, bool rDelete)
{
	mResultDir = rResultDir;
    //delete result directory
    if (rDelete)
    {
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
    else
    {
		if (boost::filesystem::exists(rResultDir))    // does p actually exist?
		{
			if (!boost::filesystem::is_directory(rResultDir))      // is p a directory?
			{
				// create result directory
				boost::filesystem::create_directory(rResultDir);
			}
		}

    }
}

//! @brief sets the minimum time step for the time integration procedure
void NuTo::TimeIntegrationBase::SetGroupNodesReactionForces(NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> rVecGroupNodesReactionForces)
{
	if (rVecGroupNodesReactionForces.GetNumColumns()!=1)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetGroupNodesReactionForces] vector (must have a single column).");
	if (rVecGroupNodesReactionForces.GetNumRows()<1)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetGroupNodesReactionForces] vector must have at least a single row.");
	mVecGroupNodesReactionForces = rVecGroupNodesReactionForces;
}

//! @brief sets the nodes, for which displacements are to be monitored
void NuTo::TimeIntegrationBase::SetOutputDispNodes(NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> rVecOutputDispNodesInt)
{
	if (rVecOutputDispNodesInt.GetNumColumns()!=1)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetOutputDispNodes] vector (must have a single column).");
	if (rVecOutputDispNodesInt.GetNumRows()<1)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetOutputDispNodes] vector must have at least a single row.");
	mVecOutputDispNodesInt = rVecOutputDispNodesInt;
}

//! @brief sets the nodes, for which displacements are to be monitored
void NuTo::TimeIntegrationBase::CalculateOutputDispNodesPtr(StructureBase& rStructure)
{
	mVecOutputDispNodesPtr.resize(mVecOutputDispNodesInt.GetNumRows());
	for (int count=0; count<mVecOutputDispNodesInt.GetNumRows(); count++)
		mVecOutputDispNodesPtr[count] = rStructure.NodeGetNodePtr(mVecOutputDispNodesInt(count,0));
}

//! @brief sets the minimum time step for the time integration procedure
void NuTo::TimeIntegrationBase::SetPlotElementGroups(NuTo::FullVector<int,Eigen::Dynamic> rPlotElementGroups)
{
	if (rPlotElementGroups.GetNumRows()<1)
		throw MechanicsException("[NuTo::TimeIntegrationBase::SetPlotElementGroups] vector must have at least a single row.");
	mPlotElementGroups = rPlotElementGroups;
}
