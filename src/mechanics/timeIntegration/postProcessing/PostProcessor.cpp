#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"

#include "boost/filesystem.hpp"

#include "base/Exception.h"

#include "mechanics/timeIntegration/postProcessing/ResultBase.h"
#include "mechanics/timeIntegration/postProcessing/ResultTime.h"
#include "mechanics/timeIntegration/postProcessing/ResultNodeDisp.h"
#include "mechanics/timeIntegration/postProcessing/ResultNodeAcceleration.h"
#include "mechanics/timeIntegration/postProcessing/ResultGroupNodeForce.h"
#include "mechanics/timeIntegration/postProcessing/ResultElementIpData.h"

using namespace NuTo;


void PostProcessor::PostProcess(const StructureOutputBlockVector& outOfBalance)
{
    if (mResultDir.length() == 0)
        throw Exception(__PRETTY_FUNCTION__, "Set the result directory first.");

    mStructure.WriteRestartFile(GetRestartFileName(), mTimeControl.GetCurrentTime());

    for (auto& result : mResults)
        result.CalculateAndAddValues(mStructure, mTimeStepResult, outOfBalance, mTimeControl.GetCurrentTime());

    if ((mTimeControl.GetCurrentTime() - mLastTimePlot) >= mMinTimeStepPlot)
    {
        for (auto& result : mResults)
            result.WriteToFile(mResultDir, mTimeStepResult);

        ExportVisualizationFiles(mTimeControl.GetCurrentTime(), mTimeStepVTK);
        mTimeStepVTK++;
        mLastTimePlot = mTimeControl.GetCurrentTime();
    }
    mTimeStepResult++;
}


void PostProcessor::AddResultNodeDisplacements(const std::string& resultName, int nodeID)
{
    mResults.push_back(new ResultNodeDisp(resultName, nodeID));
}


void PostProcessor::AddResultNodeAccelerations(const std::string& resultName, int nodeID)
{
    mResults.push_back(new ResultNodeAcceleration(resultName, nodeID));
}


void PostProcessor::AddResultTime(const std::string& resultName)
{
    mResults.push_back(new ResultTime(resultName));
}


void PostProcessor::AddResultElementIpData(const std::string& resultName, int elementID,
                                           IpData::eIpStaticDataType ipDataType)
{
    mResults.push_back(new ResultElementIpData(resultName, elementID, ipDataType));
}


void PostProcessor::AddResultGroupNodeForce(const std::string& resultName, int groupID)
{
    mResults.push_back(new ResultGroupNodeForce(resultName, groupID));
}


std::string PostProcessor::GetRestartFileName() const
{
    boost::filesystem::path restartFile(mResultDir);
    restartFile /= "LastTimeStep.restart";
    return restartFile.c_str();
}


void PostProcessor::ExportVisualizationFiles(double time, int timeStep)
{
#ifdef ENABLE_VISUALIZE
    boost::filesystem::path resultFile(mResultDir);

    if (mExportDataFileNodes == true)
    {
        resultFile /= std::string("Nodes") + std::to_string(timeStep) + std::string(".vtu");
        mStructure.ExportVtkDataFileNodes(resultFile.string());
    }

    std::stringstream timeFormatted;
    timeFormatted.width(15);
    timeFormatted.precision(12);
    timeFormatted << time;

    // plot all groups separately
    for (const int& visualizeGroup : mStructure.GetVisualizationGroups())
    {
        // plot all elements
        resultFile = mResultDir;
        resultFile /= std::string("Group") + std::to_string(visualizeGroup) + std::string("_Elements") +
                      std::to_string(timeStep) + std::string(".vtu");
        mStructure.ElementGroupExportVtkDataFile(visualizeGroup, resultFile.string());

        // write an additional pvd file
        resultFile = mResultDir;
        resultFile /= std::string("Group") + std::to_string(visualizeGroup) + std::string("_ElementsAll") +
                      std::string(".pvd");

        std::fstream file;
        if (timeStep == 0)
            file.open(resultFile.string(), std::fstream::out);
        else
            file.open(resultFile.string(), std::fstream::out | std::fstream::in | std::ios_base::ate);

        if (!file.is_open())
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Error opening file " + resultFile.string());

        std::stringstream endOfXML;
        endOfXML << "</Collection>" << std::endl;
        endOfXML << "</VTKFile>" << std::endl;
        if (timeStep == 0)
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
        file << "<DataSet timestep=\"" << timeFormatted.str() << "\" file=\"Group" << visualizeGroup
             << "_Elements" << timeStep << ".vtu\"/>" << std::endl;
        file << endOfXML.str();
        file.close();
    }

#endif // ENABLE_VISUALIZE
}


void PostProcessor::SetResultDirectory(std::string resultDir, bool deleteOld)
{
    namespace fs = boost::filesystem;

    mResultDir = resultDir;
    if (fs::exists(resultDir) and fs::is_directory(resultDir) and deleteOld)
        fs::remove_all(resultDir);
    if (not fs::exists(resultDir))
        fs::create_directory(resultDir);
    if (fs::exists(resultDir) and fs::is_regular_file(resultDir))
        throw Exception(__PRETTY_FUNCTION__, "Directory name already exists as a file.");
}
