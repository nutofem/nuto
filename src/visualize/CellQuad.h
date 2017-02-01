// $Id$

#pragma once

#include "visualize/CellBase.h"

namespace NuTo
{

//! @brief ... quadrilater cells
//! @author Stefan Eckardt, ISM
//! @date December 2009
class CellQuad: public NuTo::CellBase
{
public:
    //! @brief ... constructor
    //! @param rPoints ... point id's
    //! @param rDataTypes ... data type definitions
    CellQuad(const unsigned int *rPoints, const std::vector<VisualizeDataType>& rDataTypes);

    //! @brief ... return number of cell points
    //! @return ... number of cell points
    unsigned int GetNumPoints() const override;

    //! @brief ... return point id's
    //! @return ... array of point id's
    const unsigned int* GetPoints() const override;

    //! @brief ... returns the corresponding Vtk cell type
    //! @return ... Vtk cell type
    unsigned int GetVtkCellType() const override;
private:
    //! @brief ... point id's
    unsigned int mPoints[4];
};

}

