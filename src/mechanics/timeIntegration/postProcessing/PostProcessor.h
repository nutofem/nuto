#pragma once

#include <string>
#include <boost/ptr_container/ptr_vector.hpp>

#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/TimeControl.h"
#include "mechanics/timeIntegration/postProcessing/ResultBase.h"

namespace NuTo
{
class StructureOutputBlockVector;

//! @brief Class collecting the postprocessing options and functionality
class PostProcessor
{
public:
    PostProcessor(StructureBase& structure, const TimeControl& timeControl)
        : mStructure(structure)
        , mTimeControl(timeControl){};

    //! @brief Write out VTK file and result files
    //! @param outOfBalance Out of balance values of the independent dofs
    //! @remark `outOfBalance` here means Residual = ExternalForces - InternalForces
    void PostProcess(const StructureOutputBlockVector& outOfBalance);

    //! @brief Set the minimum time step between writing output files
    //! @param minTimeStepPlot If less than `minTimeStepPlot` time has passed between "now" and the last time an output
    //!                        file was written, postprocessing is skipped.
    void SetMinTimeStepPlot(double minTimeStepPlot)
    {
        mMinTimeStepPlot = minTimeStepPlot;
    }

    //! @brief Set the result directory
    //! @param resultDir Directory for writing results into
    //! @param deleteExisting If `deleteExisting` is set, all the content of the directory will be removed.
    void SetResultDirectory(std::string resultDir, bool deleteExisting = false);

    //! @brief Get the result directory
    std::string GetResultDirectory() const
    {
        return mResultDir;
    }

    //! @brief Get the name of the restart file
    std::string GetRestartFileName() const;

    //! @brief Monitor the accelerations of a node
    //! @param resultName Name used for the output file of the result
    //! @param nodeId ID of the node
    void AddResultNodeAccelerations(const std::string& resultName, int nodeId);

    //! @brief Monitor the displacements of a node
    //! @param resultName Name used for the output file of the result
    //! @param nodeId ID of the node
    void AddResultNodeDisplacements(const std::string& resultName, int nodeId);

    //! @brief Monitor the forces at a group of nodes
    //! @param resultName Name used for the output file of the result
    //! @param groupID ID of the group of nodes
    void AddResultGroupNodeForce(const std::string& resultName, int groupID);

    //! @brief Monitor the time
    //! @param resultName Name used for the output file of the result
    void AddResultTime(const std::string& resultName);

    //! @brief Monitor the integration point values of an element
    //! @param resultName Name used for the output file of the result
    //! @param elementID ID of the element
    //! @param ipDataType Select which IP data are to be written out
    void AddResultElementIpData(const std::string& resultName, int elementID,
                               NuTo::IpData::eIpStaticDataType ipDataType);

private:
    void ExportVisualizationFiles(double time, int timeStep);

    bool mExportDataFileNodes = true; //!< If set to true, exports a data file for the nodes

    std::string mResultDir; //!< Result directory
    boost::ptr_vector<ResultBase> mResults; //!< Specifies what to plot (displacements, reaction forces, etc.)

    int mTimeStepResult = 0; //!< Time step number is increased each time a value is added to the result matrices
    int mTimeStepVTK = 0; //!< Time step number is increased each time a VTK file is written

    double mMinTimeStepPlot = 0; //!< Minimum time between writing results
    double mLastTimePlot = -1e99; //!< Last time when a VTK file was written

    StructureBase& mStructure;

    const TimeControl& mTimeControl;
};
}
