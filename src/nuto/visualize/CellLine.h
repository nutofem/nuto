// $Id$

#ifndef CELLLINE_H_
#define CELLLINE_H_

#include "nuto/visualize/CellBase.h"

namespace NuTo
{

//! @brief ... line cells
//! @author Stefan Eckardt, ISM
//! @date 23.11.2009
class CellLine: public NuTo::CellBase
{
public:
    //! @brief ... constructor
    //! @param rPoints ... point id's
    //! @param rDataTypes ... data type definitions
    CellLine(const unsigned int *rPoints, const std::vector<VisualizeDataType>& rDataTypes);

    //! @brief ... return number of cell points
    //! @return ... number of cell points
    unsigned int GetNumPoints() const;

    //! @brief ... return point id's
    //! @return ... array of point id's
    const unsigned int* GetPoints() const;

    //! @brief ... returns the corresponding Vtk cell type
    //! @return ... Vtk cell type
    unsigned int GetVtkCellType() const;
private:
    //! @brief ... point id's
    unsigned int mPoints[2];
};

}

#endif // CELLLINE_H_ 
