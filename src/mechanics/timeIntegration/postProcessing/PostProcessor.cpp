#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"

#include "boost/filesystem.hpp"

#include "base/Timer.h"
#include "base/Exception.h"

#include "mechanics/timeIntegration/postProcessing/ResultBase.h"
#include "mechanics/timeIntegration/postProcessing/ResultTime.h"
#include "mechanics/timeIntegration/postProcessing/ResultNodeDof.h"
#include "mechanics/timeIntegration/postProcessing/ResultNodeDisp.h"
#include "mechanics/timeIntegration/postProcessing/ResultNodeAcceleration.h"
#include "mechanics/timeIntegration/postProcessing/ResultGroupNodeForce.h"
#include "mechanics/timeIntegration/postProcessing/ResultGroupNodeDof.h"
#include "mechanics/timeIntegration/postProcessing/ResultElementIpData.h"

using namespace NuTo;


void PostProcessor::PostProcess(const StructureOutputBlockVector& rOutOfBalance)
{
    if (mResultDir.length() == 0)
    {
        throw Exception(__PRETTY_FUNCTION__, "Set the result directory first.");
    }
    else
    {
        mStructure.WriteRestartFile(GetRestartFileName(), mTimeControl.GetCurrentTime());

        // perform Postprocessing
        for (auto itResult = mResultMap.begin(); itResult != mResultMap.end(); itResult++)
        {
            switch (itResult->second->GetResultType())
            {
            case eTimeIntegrationResultType::TIME:
            {
                ResultTime* resultPtr(itResult->second->AsResultTime());
                resultPtr->CalculateAndAddValues(mStructure, mTimeStepResult, mTimeControl.GetCurrentTime());
                break;
            }
            case eTimeIntegrationResultType::NODE_ACCELERATION:
            {
                ResultNodeDof* resultPtr(itResult->second->AsResultNodeDof());
                resultPtr->CalculateAndAddValues(mStructure, mTimeStepResult);
                break;
            }
            case eTimeIntegrationResultType::NODE_DISPLACEMENT:
            {
                ResultNodeDof* resultPtr(itResult->second->AsResultNodeDof());
                resultPtr->CalculateAndAddValues(mStructure, mTimeStepResult);
                break;
            }
            case eTimeIntegrationResultType::GROUP_NODE_FORCE:
            {
                ResultGroupNodeDof* resultPtr(itResult->second->AsResultGroupNodeDof());
                resultPtr->CalculateAndAddValues(mStructure, mTimeStepResult,
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
                resultPtr->CalculateAndAddValues(mStructure, mTimeStepResult);
                break;
            }
            default:
                throw Exception(__PRETTY_FUNCTION__, "Unknown component in postprocessing.");
            }
        }

        if ((mTimeControl.GetCurrentTime() - mLastTimePlot) >= mMinTimeStepPlot)
        {
            // write the results to files
            for (auto itResult = mResultMap.begin(); itResult != mResultMap.end(); itResult++)
            {
                itResult->second->WriteToFile(mResultDir, mTimeStepResult);
            }

#ifdef ENABLE_VISUALIZE
            // plot the solution vtk file
            ExportVisualizationFiles(mResultDir, mTimeControl.GetCurrentTime(), mTimeStepVTK);
#endif
            mTimeStepVTK++;
            mLastTimePlot = mTimeControl.GetCurrentTime();
        }
        mTimeStepResult++;
    }
}


int PostProcessor::AddResultNodeDisplacements(const std::string& rResultStr, int rNodeId)
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


int PostProcessor::AddResultNodeAccelerations(const std::string& rResultStr, int rNodeId)
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
    return resultNumber;
}

int PostProcessor::AddResultTime(const std::string& rResultStr)
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

int PostProcessor::AddResultElementIpData(const std::string& rResultStr, int rElementId,
                                                      IpData::eIpStaticDataType rIpDataType)
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

int PostProcessor::AddResultGroupNodeForce(const std::string& rResultStr, int rGroupNodeId)
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


std::string PostProcessor::GetRestartFileName() const
{
    boost::filesystem::path restartFile(mResultDir);
    restartFile /= "LastTimeStep.restart";
    return restartFile.c_str();
}


void PostProcessor::ExportVisualizationFiles(const std::string& rResultDir, double rTime, int rTimeStep)
{
#ifdef ENABLE_VISUALIZE
    // plot the solution vtk file
    std::stringstream ssTimeStepVTK;
    ssTimeStepVTK << rTimeStep;
    boost::filesystem::path resultFile(rResultDir);

    if (mExportDataFileNodes == true)
    {
        resultFile /= std::string("Nodes") + ssTimeStepVTK.str() + std::string(".vtu");
        mStructure.ExportVtkDataFileNodes(resultFile.string());
    }

    std::stringstream timeFormatted;
    timeFormatted.width(15);
    timeFormatted.precision(12);
    timeFormatted << rTime;

    // plot all groups separately
    for (auto const& iVisualizePair : mStructure.GetGroupVisualizeComponentsMap())
    {
        // plot all elements
        resultFile = rResultDir;
        resultFile /= std::string("Group") + std::to_string(iVisualizePair.first) + std::string("_Elements") +
                      ssTimeStepVTK.str() + std::string(".vtu");
        mStructure.ElementGroupExportVtkDataFile(iVisualizePair.first, resultFile.string());

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

void PostProcessor::SetResultDirectory(std::string rResultDir, bool rDelete)
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
