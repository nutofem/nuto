#include "boost/filesystem.hpp"

#include "base/Timer.h"

#include "mechanics/timeIntegration/TimeIntegrationBase.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/ResultElementIpData.h"
#include "mechanics/timeIntegration/ResultGroupNodeForce.h"
#include "mechanics/timeIntegration/ResultNodeDisp.h"
#include "mechanics/timeIntegration/ResultNodeAcceleration.h"
#include "mechanics/timeIntegration/ResultTime.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"

#include "mechanics/dofSubMatrixSolvers/SolverMUMPS.h"
#include "mechanics/structures/Assembler.h"

using namespace NuTo;

NuTo::TimeIntegrationBase::TimeIntegrationBase(StructureBase* rStructure)
    : mStructure(rStructure)
    , mSolver(std::make_unique<SolverMUMPS>(false))
    , mLoadVectorStatic(rStructure->GetDofStatus())
    , mLoadVectorTimeDependent(rStructure->GetDofStatus())
    , mToleranceResidual(rStructure->GetDofStatus())
    , mCallback(nullptr)
{
    ResetForNextLoad();
}

NuTo::TimeIntegrationBase::~TimeIntegrationBase()
{
}

void NuTo::TimeIntegrationBase::ResetForNextLoad()
{
    mTimeDependentLoadCase = -1;
    mTimeDependentLoadFactor.resize(0, 0);
}

const NuTo::BlockScalar& NuTo::TimeIntegrationBase::GetToleranceResidual() const
{
    return mToleranceResidual;
}

void NuTo::TimeIntegrationBase::UpdateConstraints(double rCurrentTime)
{
    mStructure->GetAssembler().ConstraintUpdateRhs(rCurrentTime);
}

void NuTo::TimeIntegrationBase::SetTimeDependentLoadCase(int rTimeDependentLoadCase,
                                                         const Eigen::MatrixXd& rTimeDependentLoadFactor)
{
    if (rTimeDependentLoadFactor.cols() != 2)
        throw MechanicsException(__PRETTY_FUNCTION__, "number of columns must be 2, first column contains the time, "
                                                      "second column contains the corresponding value.");
    if (rTimeDependentLoadFactor.rows() < 2)
        throw MechanicsException(__PRETTY_FUNCTION__, "number of rows must be at least 2.");
    if (rTimeDependentLoadFactor(0, 0) != 0)
        throw MechanicsException(__PRETTY_FUNCTION__, "the first time should always be zero.");

    // check, if the time is monotonically increasing
    for (int count = 0; count < rTimeDependentLoadFactor.rows() - 1; count++)
    {
        if (rTimeDependentLoadFactor(count, 0) >= rTimeDependentLoadFactor(count + 1, 0))
            throw MechanicsException(__PRETTY_FUNCTION__, "time has to increase monotonically.");
    }

    mTimeDependentLoadFactor = rTimeDependentLoadFactor;
    mTimeDependentLoadCase = rTimeDependentLoadCase;
}

void NuTo::TimeIntegrationBase::SetToleranceResidual(NuTo::Node::eDof rDof, double rTolerance)
{
    mToleranceResidual[rDof] = rTolerance;
}




void NuTo::TimeIntegrationBase::CalculateStaticAndTimeDependentExternalLoad()
{
    mLoadVectorStatic = StructureOutputBlockVector(mStructure->GetDofStatus(), true);
    mLoadVectorTimeDependent = StructureOutputBlockVector(mStructure->GetDofStatus(), true);

    auto tmp = mStructure->BuildGlobalExternalLoadVector();
    if (mTimeDependentLoadCase == 0)
    {
        mLoadVectorTimeDependent += tmp;
    }
    else
    {
        mLoadVectorStatic += tmp;
    }
    mStructure->GetLogger() << "Sum of loads is " << tmp.J.Export().colwise().sum() + tmp.K.Export().colwise().sum()
                            << "\n";
}




NuTo::StructureOutputBlockVector NuTo::TimeIntegrationBase::CalculateCurrentExternalLoad(double curTime)
{
    if (mTimeDependentLoadCase != -1)
    {
        if (mTimeDependentLoadFactor.rows() == 0)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "TimeDependentLoadFactor not set.");
        }
        int curStep(0);
        while (mTimeDependentLoadFactor(curStep, 0) < curTime && curStep < mTimeDependentLoadFactor.rows() - 1)
            curStep++;
        if (curStep == 0)
            curStep++;

        // extract the two data points
        double s1 = mTimeDependentLoadFactor(curStep - 1, 1);
        double s2 = mTimeDependentLoadFactor(curStep, 1);
        double t1 = mTimeDependentLoadFactor(curStep - 1, 0);
        double t2 = mTimeDependentLoadFactor(curStep, 0);

        double s = s1 + (s2 - s1) / (t2 - t1) * (curTime - t1);

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
    return mStructure->GetAssembler().GetConstraintRhs();
}

const NuTo::BlockFullVector<double>&
NuTo::TimeIntegrationBase::UpdateAndGetAndMergeConstraintRHS(double rCurrentTime, StructureOutputBlockVector& rDof_dt0)
{
    UpdateConstraints(rCurrentTime);


    rDof_dt0.K = mStructure->NodeCalculateDependentDofValues(rDof_dt0.J);
    mStructure->NodeMergeDofValues(0, rDof_dt0);

    mStructure->ElementTotalUpdateTmpStaticData();
    return mStructure->GetAssembler().GetConstraintRhs();
}

int NuTo::TimeIntegrationBase::AddResultNodeDisplacements(const std::string& rResultStr, int rNodeId)
{
    // find unused integer id
    int resultNumber(mResultMap.size());
    boost::ptr_map<int, ResultBase>::iterator it = mResultMap.find(resultNumber);
    while (it != mResultMap.end())
    {
        resultNumber++;
        it = mResultMap.find(resultNumber);
    }

    mResultMap.insert(resultNumber, new ResultNodeDisp(rResultStr, rNodeId));

    return resultNumber;
}

int NuTo::TimeIntegrationBase::AddResultNodeAccelerations(const std::string& rResultStr, int rNodeId)
{
    // find unused integer id
    int resultNumber(mResultMap.size());
    boost::ptr_map<int, ResultBase>::iterator it = mResultMap.find(resultNumber);
    while (it != mResultMap.end())
    {
        resultNumber++;
        it = mResultMap.find(resultNumber);
    }

    mResultMap.insert(resultNumber, new ResultNodeAcceleration(rResultStr, rNodeId));
    mMergeActiveDofValuesOrder2 = true;
    return resultNumber;
}

int NuTo::TimeIntegrationBase::AddResultTime(const std::string& rResultStr)
{
    // find unused integer id
    int resultNumber(mResultMap.size());
    boost::ptr_map<int, ResultBase>::iterator it = mResultMap.find(resultNumber);
    while (it != mResultMap.end())
    {
        resultNumber++;
        it = mResultMap.find(resultNumber);
    }

    mResultMap.insert(resultNumber, new ResultTime(rResultStr));

    return resultNumber;
}

int NuTo::TimeIntegrationBase::AddResultElementIpData(const std::string& rResultStr, int rElementId,
                                                      NuTo::IpData::eIpStaticDataType rIpDataType)
{
    // find unused integer id
    int resultNumber(mResultMap.size());
    boost::ptr_map<int, ResultBase>::iterator it = mResultMap.find(resultNumber);
    while (it != mResultMap.end())
    {
        resultNumber++;
        it = mResultMap.find(resultNumber);
    }

    mResultMap.insert(resultNumber, new ResultElementIpData(rResultStr, rElementId, rIpDataType));

    return resultNumber;
}

int NuTo::TimeIntegrationBase::AddResultGroupNodeForce(const std::string& rResultStr, int rGroupNodeId)
{
    // find unused integer id
    int resultNumber(mResultMap.size());
    boost::ptr_map<int, ResultBase>::iterator it = mResultMap.find(resultNumber);
    while (it != mResultMap.end())
    {
        resultNumber++;
        it = mResultMap.find(resultNumber);
    }

    mResultMap.insert(resultNumber, new ResultGroupNodeForce(rResultStr, rGroupNodeId));

    return resultNumber;
}

std::array<StructureOutputBlockVector, 3> NuTo::TimeIntegrationBase::ExtractDofValues() const
{
    const auto& dofStatus = mStructure->GetDofStatus();
    std::array<StructureOutputBlockVector, 3> dofValues = {StructureOutputBlockVector(dofStatus, true),
                                                           StructureOutputBlockVector(dofStatus, true),
                                                           StructureOutputBlockVector(dofStatus, true)};
    dofValues[0] = mStructure->NodeExtractDofValues(0);

    if (mStructure->GetNumTimeDerivatives() >= 1)
        dofValues[1] = mStructure->NodeExtractDofValues(1);

    if (mStructure->GetNumTimeDerivatives() >= 2)
        dofValues[2] = mStructure->NodeExtractDofValues(2);
    return dofValues;
}

double NuTo::TimeIntegrationBase::CalculateNorm(const BlockFullVector<double>& rResidual) const
{
    double norm = 0;
    for (auto rDofType : rResidual.GetDofStatus().GetActiveDofTypes())
        norm += rResidual[rDofType].norm();

    return norm;
}

std::string NuTo::TimeIntegrationBase::GetRestartFileName() const
{
    boost::filesystem::path restartFile(mResultDir);
    restartFile /= "LastTimeStep.restart";
    return restartFile.c_str();
}

void NuTo::TimeIntegrationBase::PostProcess(const StructureOutputBlockVector& rOutOfBalance)
{
    Timer timer(__FUNCTION__, GetShowTime(), mStructure->GetLogger());

    if (mResultDir.length() == 0)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Set the result directory first.");
    }
    else
    {
        mStructure->WriteRestartFile(GetRestartFileName(), mTime);

        // perform Postprocessing
        for (auto itResult = mResultMap.begin(); itResult != mResultMap.end(); itResult++)
        {
            switch (itResult->second->GetResultType())
            {
            case eTimeIntegrationResultType::TIME:
            {
                ResultTime* resultPtr(itResult->second->AsResultTime());
                resultPtr->CalculateAndAddValues(*mStructure, mTimeStepResult, mTime);
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
                resultPtr->CalculateAndAddValues(*mStructure, mTimeStepResult,
                                                 rOutOfBalance.J[Node::eDof::DISPLACEMENTS],
                                                 rOutOfBalance.K[Node::eDof::DISPLACEMENTS]);
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

        if ((mTime - mLastTimePlot) >= mMinTimeStepPlot)
        {
            // write the results to files
            for (auto itResult = mResultMap.begin(); itResult != mResultMap.end(); itResult++)
            {
                itResult->second->WriteToFile(mResultDir, mTimeStepResult);
            }

#ifdef ENABLE_VISUALIZE
            // plot the solution vtk file
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

void NuTo::TimeIntegrationBase::SetActiveDofsCalculationStep(int rStepNum,
                                                             const std::set<NuTo::Node::eDof>& rActiveDofs)
{
    mStepActiveDofs[rStepNum] = rActiveDofs;
}

void NuTo::TimeIntegrationBase::Info() const
{
}

void NuTo::TimeIntegrationBase::SetResultDirectory(std::string rResultDir, bool rDelete)
{
    mResultDir = rResultDir;
    // delete result directory
    if (rDelete)
    {
        if (boost::filesystem::exists(rResultDir)) // does p actually exist?
        {
            if (boost::filesystem::is_directory(rResultDir)) // is p a directory?
            {
                boost::filesystem::remove_all(rResultDir);
            }
        }
        // create result directory
        boost::filesystem::create_directory(rResultDir);
    }
    else
    {
        if (boost::filesystem::exists(rResultDir)) // does p actually exist?
        {
            if (!boost::filesystem::is_directory(rResultDir)) // is p a directory?
            {
                // create result directory
                boost::filesystem::create_directory(rResultDir);
            }
        }
    }
}


void NuTo::TimeIntegrationBase::ExportVisualizationFiles(const std::string& rResultDir, double rTime, int rTimeStep)
{
#ifdef ENABLE_VISUALIZE
    // plot the solution vtk file
    std::stringstream ssTimeStepVTK;
    ssTimeStepVTK << rTimeStep;
    boost::filesystem::path resultFile(rResultDir);

    if (mExportDataFileNodes == true)
    {
        resultFile /= std::string("Nodes") + ssTimeStepVTK.str() + std::string(".vtu");
        mStructure->ExportVtkDataFileNodes(resultFile.string());
    }

    std::stringstream timeFormatted;
    timeFormatted.width(15);
    timeFormatted.precision(12);
    timeFormatted << rTime;

    // plot all groups separately
    for (auto const& iVisualizePair : mStructure->GetGroupVisualizeComponentsMap())
    {
        // plot all elements
        resultFile = rResultDir;
        resultFile /= std::string("Group") + std::to_string(iVisualizePair.first) + std::string("_Elements") +
                      ssTimeStepVTK.str() + std::string(".vtu");
        mStructure->ElementGroupExportVtkDataFile(iVisualizePair.first, resultFile.string());

        // write an additional pvd file
        resultFile = rResultDir;
        resultFile /= std::string("Group") + std::to_string(iVisualizePair.first) + std::string("_ElementsAll") +
                      std::string(".pvd");

        std::fstream file;
        if (rTimeStep == 0)
        {
            file.open(resultFile.string(), std::fstream::out);
        }
        else
        {
            file.open(resultFile.string(), std::fstream::out | std::fstream::in | std::ios_base::ate);
        }
        if (!file.is_open())
        {
            throw NuTo::MechanicsException(
                    std::string("[NuTo::TimeIntegrationBase::ExportVisualizationFiles] Error opening file ") +
                    resultFile.string());
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
        }
        else
        {
            // delete the last part of the xml file
            file.seekp(-endOfXML.str().length(), std::ios_base::end);
        }
        file << "<DataSet timestep=\"" << timeFormatted.str() << "\" file=\"Group" << iVisualizePair.first
             << "_Elements" << rTimeStep << ".vtu\"/>" << std::endl;
        file << endOfXML.str();
        file.close();
    }

#endif // ENABLE_VISUALIZE
}

bool TimeIntegrationBase::GetShowTime() const
{
    return mShowTime;
}

void TimeIntegrationBase::SetShowTime(bool showTime)
{
    mShowTime = showTime;
}
