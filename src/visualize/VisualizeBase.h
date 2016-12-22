// $Id$

#pragma once

#include <string>

namespace NuTo
{

//! @brief ... base class for visualize objects
//! @author Stefan Eckardt, ISM
//! @date November 2009
class VisualizeBase
{
public:

    //! @brief ... export to Vtk datafile
    //! @param rFilename ... filename
    void ExportVtkDataFile(const std::string& rFilename) const;
};

}

