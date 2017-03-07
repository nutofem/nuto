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


#include "base/Timer.h"

#include "mechanics/timeIntegration/TimeIntegrationBase.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/ResultElementIpData.h"
#include "mechanics/timeIntegration/ResultGroupNodeForce.h"
#include "mechanics/timeIntegration/ResultNodeDisp.h"
#include "mechanics/timeIntegration/ResultNodeAcceleration.h"
#include "mechanics/timeIntegration/ResultTime.h"
#include "mechanics/timeIntegration/TimeDependencyFunction.h"
#include "mechanics/timeIntegration/TimeDependencyMatrix.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"

#include "mechanics/dofSubMatrixSolvers/SolverMUMPS.h"



NuTo::TimeIntegrationBase::TimeIntegrationBase(StructureBase* rStructure) :
        NuTo::NuToObject::NuToObject(),
        mStructure(rStructure),
        mSolver(std::make_unique<SolverMUMPS>(false)),
        mTimeDependentConstraint(-1),
        mTimeDependentLoadCase(-1),
        mLoadVectorStatic(rStructure->GetDofStatus()),
        mLoadVectorTimeDependent(rStructure->GetDofStatus()),
        mTime(0.),
        mAutomaticTimeStepping(false),
        mTimeStep(0),
        mMaxTimeStep(1),
        mMinTimeStep(0),
        mMergeActiveDofValuesOrder1(true),
        mMergeActiveDofValuesOrder2(false),
        mCheckCoefficientMatrix(false),
        mToleranceResidual(rStructure->GetDofStatus()),
        mLoadStep(0),
        mTimeStepResult(0),
        mTimeStepVTK(0),
        mMinTimeStepPlot(0),
        mLastTimePlot(-1e99),
        mIterationCount(0),
        mCallback(nullptr)
{
    ResetForNextLoad();
}

NuTo::TimeIntegrationBase::~TimeIntegrationBase()
{}

//! @brief sets the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
void NuTo::TimeIntegrationBase::ResetForNextLoad()
{
    mTimeDependentConstraint = -1;
    mTimeDependentConstraintFactor.resize(0,0);
    mTimeDependentLoadCase = -1;
    mTimeDependentLoadFactor.resize(0,0);
}


//! @brief Adds the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
//! @param rTimeDependentConstraint ... constraint, whose rhs is increased as a function of time
//! @param rTimeDependentConstraintFactor ... first row time, rhs of the constraint (linear interpolation in between afterwards linear extrapolation)
void NuTo::TimeIntegrationBase::AddTimeDependentConstraint(int rTimeDependentConstraint, const Eigen::MatrixXd &rTimeDependentConstraintFactor)
{
    if (rTimeDependentConstraintFactor.cols()!=2)
        throw MechanicsException(__PRETTY_FUNCTION__, "number of columns must be 2, first column contains the time, second column contains the corresponding rhs.");
    if (rTimeDependentConstraintFactor.rows()<2)
        throw MechanicsException(__PRETTY_FUNCTION__, "number of rows must be at least 2.");
    if (rTimeDependentConstraintFactor(0,0)!=0)
        throw MechanicsException(__PRETTY_FUNCTION__, "the first time should always be zero.");
    //check, if the time is monotonically increasing
    for (int count=0; count<rTimeDependentConstraintFactor.rows()-1; count++)
    {
        if (rTimeDependentConstraintFactor(count,0)>=rTimeDependentConstraintFactor(count+1,0))
            throw MechanicsException(__PRETTY_FUNCTION__, "time has to increase monotonically.");
    }

    mMapTimeDependentConstraint.insert(std::pair<int,std::shared_ptr<TimeDependencyBase>> (rTimeDependentConstraint, std::make_shared<TimeDependencyMatrix>(rTimeDependentConstraintFactor)));
}


//! @brief Adds the delta rhs of the constrain equation whose RHS is incrementally increased in each load step / time step
//! @param rTimeDependentConstraint ... constraint, whose rhs is increased as a function of time
//! @param rTimeDependentConstraintFunction ... function that calculates the time dependent constraint factor for the current time step
void NuTo::TimeIntegrationBase::AddTimeDependentConstraintFunction(int rTimeDependentConstraint, const std::function<double (double rTime)>& rTimeDependentConstraintFunction)
{
    mMapTimeDependentConstraint.insert(std::pair<int,std::shared_ptr<TimeDependencyBase>> (rTimeDependentConstraint, std::make_shared<TimeDependencyFunction>(rTimeDependentConstraintFunction)));
}

const NuTo::BlockScalar& NuTo::TimeIntegrationBase::GetToleranceResidual() const
{
    return mToleranceResidual;
}

//! @brief Updates the Rhs for all constraints
//! @param rCurrentTime ... current time
//! @remark remove the second argument rDof
void NuTo::TimeIntegrationBase::UpdateConstraints(double rCurrentTime)
{
    for(auto itTDC : mMapTimeDependentConstraint)
        mStructure->ConstraintSetRHS(itTDC.first, itTDC.second->GetTimeDependentFactor(rCurrentTime));
}


//! @brief sets a scalar time dependent multiplication factor for the external loads
//! @param rTimeDependentLoadFactor ... first row time, second row scalar factor to calculate the external load (linear interpolation in between, afterwards constant)
void NuTo::TimeIntegrationBase::SetTimeDependentLoadCase(int rTimeDependentLoadCase, const Eigen::MatrixXd& rTimeDependentLoadFactor)
{
    if (rTimeDependentLoadFactor.cols()!=2)
        throw MechanicsException(__PRETTY_FUNCTION__, "number of columns must be 2, first column contains the time, second column contains the corresponding value.");
    if (rTimeDependentLoadFactor.rows()<2)
        throw MechanicsException(__PRETTY_FUNCTION__, "number of rows must be at least 2.");
    if (rTimeDependentLoadFactor(0,0)!=0)
        throw MechanicsException(__PRETTY_FUNCTION__, "the first time should always be zero.");

    //check, if the time is monotonically increasing
    for (int count=0; count<rTimeDependentLoadFactor.rows()-1; count++)
    {
        if (rTimeDependentLoadFactor(count,0)>=rTimeDependentLoadFactor(count+1,0))
            throw MechanicsException(__PRETTY_FUNCTION__, "time has to increase monotonically.");
    }

    mTimeDependentLoadFactor = rTimeDependentLoadFactor;
    mTimeDependentLoadCase = rTimeDependentLoadCase;
}

//! @brief apply the new rhs of the constraints as a function of the current time delta
double NuTo::TimeIntegrationBase::CalculateTimeDependentConstraintFactor(double curTime)
{
    //calculate the two corresponding time steps between which a linear interpolation is performed
    if (mTimeDependentConstraintFactor.rows()!=0)
    {
        int curStep(0);
        while (mTimeDependentConstraintFactor(curStep,0)<curTime && curStep<mTimeDependentConstraintFactor.rows()-1)
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




//! @brief Sets the residual tolerance for a specific DOF
//! param rDof: degree of freedom
//! param rTolerance: tolerance
void NuTo::TimeIntegrationBase::SetToleranceResidual(NuTo::Node::eDof rDof, double rTolerance)
{
    mToleranceResidual[rDof] = rTolerance;
}





//! @brief calculate the external force vector (mStatic and m) as a function of time delta
void NuTo::TimeIntegrationBase::CalculateStaticAndTimeDependentExternalLoad()
{
    mLoadVectorStatic        = StructureOutputBlockVector(mStructure->GetDofStatus(), true);
    mLoadVectorTimeDependent = StructureOutputBlockVector(mStructure->GetDofStatus(), true);

    for (int iLoadCase = 0; iLoadCase< mStructure->GetNumLoadCases(); ++iLoadCase)
    {
        auto tmp = mStructure->BuildGlobalExternalLoadVector(iLoadCase);
        mStructure->GetLogger()<<"TIB_CEL1, mTimeDeoendentLoadCase = " << mTimeDependentLoadCase << "\n";
        if (iLoadCase == mTimeDependentLoadCase)
        {
            mStructure->GetLogger() << "TIB TimeDependent \n";
            mLoadVectorTimeDependent += tmp;
        }
        else
        {
            mStructure->GetLogger() << "TIB Static \n";
            mLoadVectorStatic += tmp;
        }
        mStructure->GetLogger() << "sum of loads for loadcase " << iLoadCase << " is " << 
            tmp.J.Export().colwise().sum() + tmp.K.Export().colwise().sum() << "\n";
    }
}

//! @brief calculate the current external force as a function of time delta
//! @param curTime ... current time in the load step
//! @return ... external load vector
NuTo::StructureOutputBlockVector NuTo::TimeIntegrationBase::CalculateCurrentExternalLoad(double curTime)
{
    if (mTimeDependentLoadCase!=-1)
    {
        if (mTimeDependentLoadFactor.rows()==0)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "TimeDependentLoadFactor not set.");
        }
        int curStep(0);
        while (mTimeDependentLoadFactor(curStep,0)<curTime && curStep<mTimeDependentLoadFactor.rows()-1)
            curStep++;
        if (curStep==0)
            curStep++;

        //extract the two data points
        double s1 = mTimeDependentLoadFactor(curStep-1,1);
        double s2 = mTimeDependentLoadFactor(curStep,1);
        double t1 = mTimeDependentLoadFactor(curStep-1,0);
        double t2 = mTimeDependentLoadFactor(curStep,0);

        double s =  s1 + (s2-s1)/(t2-t1) * (curTime-t1);

        return mLoadVectorStatic + mLoadVectorTimeDependent * s;
    }
    else
    {
        return mLoadVectorStatic;
    }
}

const NuTo::BlockFullVector<double>& NuTo::TimeIntegrationBase::UpdateAndGetConstraintRHS(double rCurrentTime)
{
    UpdateConstraints(rCurrentTime);
    return mStructure->ConstraintGetRHSAfterGaussElimination();
}

const NuTo::BlockFullVector<double>& NuTo::TimeIntegrationBase::UpdateAndGetAndMergeConstraintRHS(double rCurrentTime, StructureOutputBlockVector& rDof_dt0)
{
    UpdateConstraints(rCurrentTime);


    rDof_dt0.K = mStructure->NodeCalculateDependentDofValues(rDof_dt0.J);
    mStructure->NodeMergeDofValues(0, rDof_dt0);



    mStructure->ElementTotalUpdateTmpStaticData();

    return mStructure->ConstraintGetRHSAfterGaussElimination();
}

//! @brief monitor the displacements of a node
//! @param rNodeId id of the node
//! @param rResultId string identifying the result, this is used for the output file
//! @return id of the result, so that it could be modified afterwards
int NuTo::TimeIntegrationBase::AddResultNodeDisplacements(const std::string& rResultStr, int rNodeId )
{
    //find unused integer id
    int resultNumber(mResultMap.size());
    boost::ptr_map<int,ResultBase>::iterator it = mResultMap.find(resultNumber);
    while (it!=mResultMap.end())
    {
        resultNumber++;
        it = mResultMap.find(resultNumber);
    }

    mResultMap.insert(resultNumber, new ResultNodeDisp(rResultStr,rNodeId));

    return resultNumber;
}

//! @brief monitor the accelerations of a node
//! @param rNodeId id of the node
//! @param rResultId string identifying the result, this is used for the output file
//! @return id of the result, so that it could be modified afterwards
int NuTo::TimeIntegrationBase::AddResultNodeAccelerations(const std::string& rResultStr, int rNodeId )
{
    //find unused integer id
    int resultNumber(mResultMap.size());
    boost::ptr_map<int,ResultBase>::iterator it = mResultMap.find(resultNumber);
    while (it!=mResultMap.end())
    {
        resultNumber++;
        it = mResultMap.find(resultNumber);
    }

    mResultMap.insert(resultNumber, new ResultNodeAcceleration(rResultStr,rNodeId));
    mMergeActiveDofValuesOrder2 = true;
    return resultNumber;
}

//! @brief monitor the time
//! @param rResultId string identifying the result, this is used for the output file
//! @return id of the result, so that it could be modified afterwards
int NuTo::TimeIntegrationBase::AddResultTime(const std::string& rResultStr)
{
    //find unused integer id
    int resultNumber(mResultMap.size());
    boost::ptr_map<int,ResultBase>::iterator it = mResultMap.find(resultNumber);
    while (it!=mResultMap.end())
    {
        resultNumber++;
        it = mResultMap.find(resultNumber);
    }

    mResultMap.insert(resultNumber, new ResultTime(rResultStr));

    return resultNumber;
}


int NuTo::TimeIntegrationBase::AddResultElementIpData(const std::string& rResultStr, int rElementId, NuTo::IpData::eIpStaticDataType rIpDataType)
{
    //find unused integer id
    int resultNumber(mResultMap.size());
    boost::ptr_map<int,ResultBase>::iterator it = mResultMap.find(resultNumber);
    while (it!=mResultMap.end())
    {
        resultNumber++;
        it = mResultMap.find(resultNumber);
    }

    mResultMap.insert(resultNumber, new ResultElementIpData(rResultStr, rElementId, rIpDataType));

    return resultNumber;
}


//! @brief monitor the time
//! @param rResultId string identifying the result, this is used for the output file
//! @param rGroupNodeId group id of the node group, for which the reaction forces (out of balance forces) should be calculated
//! @return id of the result, so that it could be modified afterwards
int NuTo::TimeIntegrationBase::AddResultGroupNodeForce(const std::string& rResultStr,int rGroupNodeId)
{
    //find unused integer id
    int resultNumber(mResultMap.size());
    boost::ptr_map<int,ResultBase>::iterator it = mResultMap.find(resultNumber);
    while (it!=mResultMap.end())
    {
        resultNumber++;
        it = mResultMap.find(resultNumber);
    }

    mResultMap.insert(resultNumber, new ResultGroupNodeForce(rResultStr, rGroupNodeId));

    return resultNumber;
}

//! @brief extracts all dof values
//! @param rDof_dt0 ... 0th time derivative
//! @param rDof_dt1 ... 1st time derivative
//! @param rDof_dt2 ... 2nd time derivative
void NuTo::TimeIntegrationBase::ExtractDofValues(StructureOutputBlockVector& rDof_dt0, StructureOutputBlockVector& rDof_dt1, StructureOutputBlockVector& rDof_dt2) const
{
    rDof_dt0 = mStructure->NodeExtractDofValues(0);

    if (mStructure->GetNumTimeDerivatives() >= 1)
        rDof_dt1 = mStructure->NodeExtractDofValues(1);

    if (mStructure->GetNumTimeDerivatives() >= 2)
        rDof_dt2 = mStructure->NodeExtractDofValues(2);

}

//! @brief calculates the norm of the residual, can include weighting
//! @param rResidual ... residual
double NuTo::TimeIntegrationBase::CalculateNorm(const BlockFullVector<double>& rResidual) const
{
    double norm = 0;
    for (auto rDofType : rResidual.GetDofStatus().GetActiveDofTypes())
        norm += rResidual[rDofType].norm();

    return norm;
}

//! @brief postprocess (nodal dofs etc. and visualize a vtk file)
//! @param rOutOfBalance ... out of balance values of the independent dofs (for disp dofs, this is the out of balance force)
void NuTo::TimeIntegrationBase::PostProcess(const StructureOutputBlockVector& rOutOfBalance)
{
    Timer timer(__FUNCTION__, GetShowTime(), mStructure->GetLogger());

    if (mResultDir.length()==0)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Set the result directory first.");
    }
    else
    {
        //perform Postprocessing
        for (auto itResult=mResultMap.begin(); itResult!=mResultMap.end(); itResult++)
        {
            switch (itResult->second->GetResultType())
            {
            case eTimeIntegrationResultType::TIME:
            {
                ResultTime* resultPtr(itResult->second->AsResultTime());
                resultPtr->CalculateAndAddValues(*mStructure, mTimeStepResult,mTime);
                break;
            }
            case eTimeIntegrationResultType::NODE_ACCELERATION:
            {
                ResultNodeDof* resultPtr(itResult->second->AsResultNodeDof());
                resultPtr->CalculateAndAddValues(*mStructure, mTimeStepResult);
                break;
            }
            case eTimeIntegrationResultType::NODE_DISPLACEMENT:
            {
                ResultNodeDof* resultPtr(itResult->second->AsResultNodeDof());
                resultPtr->CalculateAndAddValues(*mStructure, mTimeStepResult);
                break;
            }
            case eTimeIntegrationResultType::GROUP_NODE_FORCE:
            {
                ResultGroupNodeDof* resultPtr(itResult->second->AsResultGroupNodeDof());
                resultPtr->CalculateAndAddValues(*mStructure, mTimeStepResult, rOutOfBalance.J[Node::eDof::DISPLACEMENTS],rOutOfBalance.K[Node::eDof::DISPLACEMENTS]);
                break;

            }
            case eTimeIntegrationResultType::ELEMENT_IP_STRESS:
            case eTimeIntegrationResultType::ELEMENT_IP_STRAIN:
            case eTimeIntegrationResultType::ELEMENT_IP_DAMAGE:
            case eTimeIntegrationResultType::ELEMENT_IP_BOND_STRESS:
            case eTimeIntegrationResultType::ELEMENT_IP_SLIP:
            {
                ResultElementIpData* resultPtr(itResult->second->AsResultElementIpData());
                resultPtr->CalculateAndAddValues(*mStructure, mTimeStepResult);
                break;
            }
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Unknown component in postprocessing.");
            }
        }

        if ((mTime-mLastTimePlot)>=mMinTimeStepPlot)
        {
            //write the results to files
            for (auto itResult=mResultMap.begin(); itResult!=mResultMap.end(); itResult++)
            {
                itResult->second->WriteToFile(mResultDir,mTimeStepResult);
            }

#ifdef ENABLE_VISUALIZE
            //plot the solution vtk file
            ExportVisualizationFiles(mResultDir, mTime, mTimeStepVTK);
#endif
            mTimeStepVTK++;
            mLastTimePlot = mTime;
        }
        mTimeStepResult++;
    }
}

void NuTo::TimeIntegrationBase::AddCalculationStep(const std::set<NuTo::Node::eDof>& rActiveDofs)
{
    mStepActiveDofs.push_back(rActiveDofs);
}


void NuTo::TimeIntegrationBase::SetNumCalculationSteps(int rNumSteps)
{
    mStepActiveDofs.resize(rNumSteps);
}

void NuTo::TimeIntegrationBase::SetActiveDofsCalculationStep(int rStepNum, const std::set<NuTo::Node::eDof>& rActiveDofs)
{
    mStepActiveDofs[rStepNum] = rActiveDofs;
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
    & BOOST_SERIALIZATION_NVP(mTimeDependentConstraintFactor)
    //       & BOOST_SERIALIZATION_NVP(mMapTimeDependentConstraint)
    & BOOST_SERIALIZATION_NVP(mTimeDependentLoadCase)
    & BOOST_SERIALIZATION_NVP(mTimeDependentLoadFactor)
    //       & BOOST_SERIALIZATION_NVP(mLoadVectorStatic)
    //       & BOOST_SERIALIZATION_NVP(mLoadVectorTimeDependent)
    & BOOST_SERIALIZATION_NVP(mTime)
    & BOOST_SERIALIZATION_NVP(mTimeStepResult)
    & BOOST_SERIALIZATION_NVP(mTimeStepVTK)
    & BOOST_SERIALIZATION_NVP(mLoadStep)
    & BOOST_SERIALIZATION_NVP(mTimeStep)
    & BOOST_SERIALIZATION_NVP(mMaxTimeStep)
    & BOOST_SERIALIZATION_NVP(mMinTimeStep)
    & BOOST_SERIALIZATION_NVP(mLoadStep)
    & BOOST_SERIALIZATION_NVP(mMinTimeStepPlot)
    & BOOST_SERIALIZATION_NVP(mLastTimePlot)
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
void NuTo::TimeIntegrationBase::SetPlotElementGroups(std::vector<int> rPlotElementGroups)
{
    if (rPlotElementGroups.empty())
        throw MechanicsException(__PRETTY_FUNCTION__, "vector must have at least a single row.");
    mPlotElementGroups = rPlotElementGroups;
}

void NuTo::TimeIntegrationBase::ExportVisualizationFiles(const std::string& rResultDir, double rTime, int rTimeStep)
{
#ifdef ENABLE_VISUALIZE
    //plot the solution vtk file
    std::stringstream ssTimeStepVTK;
    ssTimeStepVTK << rTimeStep;
    boost::filesystem::path resultFile(rResultDir);

    if (mExportDataFileNodes==true)
    {
        resultFile /= std::string("Nodes") + ssTimeStepVTK.str() + std::string(".vtu");
        mStructure->ExportVtkDataFileNodes(resultFile.string(), true);
    }

    std::stringstream timeFormatted;
    timeFormatted.width(15);
    timeFormatted.precision(12);
    timeFormatted << rTime;

    //plot all groups separately
    for (auto const & iVisualizePair : mStructure->GetGroupVisualizeComponentsMap())
    {
        //plot all elements
        resultFile = rResultDir;
        resultFile /= std::string("Group") + std::to_string(iVisualizePair.first) + std::string("_Elements") + ssTimeStepVTK.str() + std::string(".vtu");
        mStructure->ElementGroupExportVtkDataFile(iVisualizePair.first, resultFile.string(), true);

        //write an additional pvd file
        resultFile = rResultDir;
        resultFile /= std::string("Group") + std::to_string(iVisualizePair.first) + std::string("_ElementsAll") + std::string(".pvd");

        std::fstream file;
        if (rTimeStep == 0)
        {
            file.open(resultFile.string(), std::fstream::out);
        } else
        {
            file.open(resultFile.string(), std::fstream::out | std::fstream::in | std::ios_base::ate);
        }
        if (!file.is_open())
        {
            throw NuTo::MechanicsException(std::string("[NuTo::TimeIntegrationBase::ExportVisualizationFiles] Error opening file ") + resultFile.string());
        }
        std::stringstream endOfXML;
        endOfXML << "</Collection>" << std::endl;
        endOfXML << "</VTKFile>" << std::endl;
        if (rTimeStep == 0)
        {
            // header /////////////////////////////////////////////////////////////////
            file << "<?xml version=\"1.0\"?>" << std::endl;
            file << "<VTKFile type=\"Collection\">" << std::endl;
            file << "<Collection>" << std::endl;
        } else
        {
            //delete the last part of the xml file
            file.seekp(-endOfXML.str().length(), std::ios_base::end);
        }
        file << "<DataSet timestep=\"" << timeFormatted.str() << "\" file=\"Group" << iVisualizePair.first << "_Elements" << rTimeStep << ".vtu\"/>" << std::endl;
        file << endOfXML.str();
        file.close();
    }

#endif //ENABLE_VISUALIZE
}
