// $Id$

#ifndef VISUALIZEBASE_H_
#define VISUALIZEBASE_H_

#include <string>

namespace NuTo
{

//! @brief ... base class for visualize objects
//! @author Stefan Eckardt, ISM
//! @date November 2009
class VisualizeBase
{
public:
    enum eVisualizeWhat
    {
        DISPLACEMENTS,            //!< visualize displacements
        ENGINEERING_STRESS,       //!< visualize engineering stress tensor
        ENGINEERING_STRAIN        //!< visualize engineering strain tensor
    };

    //! @brief ... export to Vtk datafile
    //! @param rFilename ... filename
    void ExportVtkDataFile(const std::string& rFilename) const;
};

}

#endif // VISUALIZEBASE_H_ 
