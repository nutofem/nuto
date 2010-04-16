// $Id$

#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeException.h"

// export to Vtk
void NuTo::VisualizeBase::ExportVtkDataFile(const std::string& rFilename) const
{
    throw NuTo::VisualizeException("[NuTo::VisualizeBase::ExportVtkDatafile] Export to Vtk data file is not supported");
}

