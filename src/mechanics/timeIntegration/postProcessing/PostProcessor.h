#pragma once

#include <string>
#include <boost/ptr_container/ptr_vector.hpp>

#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/TimeControl.h"
#include "mechanics/timeIntegration/postProcessing/ResultBase.h"

namespace NuTo
{
class StructureOutputBlockVector;

class PostProcessor
{
public:
    PostProcessor(StructureBase& structure, const TimeControl& timeControl)
        : mStructure(structure)
        , mTimeControl(timeControl){};

    //! @brief postprocess (nodal dofs etc. and visualize a vtk file)
    //! @param rOutOfBalance ... out of balance values of the independent dofs (for disp dofs, this is the out of
    //! balance force)
    //! @remark rOutOfBalance here means Residual = ExternalForces - InternalForces
    void PostProcess(const StructureOutputBlockVector& rOutOfBalance);

    //! @brief sets the minimum time step for the time integration procedure
    void SetMinTimeStepPlot(double rMinTimeStepPlot)
    {
        mMinTimeStepPlot = rMinTimeStepPlot;
    }

    //! @brief sets the result directory
    //! @param if delete is set, all the content of the directory will be removed
    void SetResultDirectory(std::string rResultDir, bool rDelete = false);

    //! @brief getter for the result directory
    std::string GetResultDirectory() const
    {
        return mResultDir;
    }

    //! @brief returns the name of the restart file
    std::string GetRestartFileName() const;

    //! @brief monitor the accelerations of a node
    //! @param rNodeId id of the node
    //! @param rResultId string identifying the result, this is used for the output file
    void AddResultNodeAccelerations(const std::string& rResultStr, int rNodeId);

    //! @brief monitor the displacements of a node
    //! @param rNodeId id of the node
    //! @param rResultId string identifying the result, this is used for the output file
    void AddResultNodeDisplacements(const std::string& rResultStr, int rNodeId);

    //! @brief monitor the time
    //! @param rResultId string identifying the result, this is used for the output file
    //! @param rGroupNodeId group id of the node group, for which the reaction forces (out of balance forces) should be
    //! calculated
    void AddResultGroupNodeForce(const std::string& rResultStr, int rGroupNodeId);

    //! @brief monitor the time
    //! @param rResultId string identifying the result, this is used for the output file
    void AddResultTime(const std::string& rResultStr);

    //! @brief monitor the integration point values in an element
    //! @param rResultId string identifying the result, this is used for the output file
    //! @param rElementId id of the element to be monitored
    void AddResultElementIpData(const std::string& rResultStr, int rElementId,
                               NuTo::IpData::eIpStaticDataType rIpDataType);

private:
    void ExportVisualizationFiles(const std::string& rResultDir, double rTime, int timeStep);

    bool mExportDataFileNodes = true; //!< if set to true, exports a data file for the nodes

    std::string mResultDir; //!< result directory
    boost::ptr_vector<ResultBase> mResults; //!< specifies what to plot (displacements, reaction forces, etc.)

    int mTimeStepResult = 0; //!< time step number is increased each time a value is added to the result matrices
    int mTimeStepVTK = 0; //!< time step number is increased each time a vtk file is extracted
    double mMinTimeStepPlot = 0; //!< if the time between the current time step and the previous plotted step is larger
    //! than mMaxDeltaTimeStepPlot a vtk plot is performed
    double mLastTimePlot = -1e99; //!< last time when a vtk file was plotted

    StructureBase& mStructure;

    const TimeControl& mTimeControl;
};
}
