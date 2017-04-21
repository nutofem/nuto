// $Id$
#pragma once

#include <vector>
#include <eigen3/Eigen/Core>


namespace NuTo
{

//! @brief ... base class for visualization cells
//! @author Stefan Eckardt, ISM
//! @date November 2009
class CellBase
{
public:

    //! @brief constructor
    //! @param numData ... number of different data fields
    CellBase(int numData);

    //! @brief ... return number of cell points
    //! @return ... number of cell points
    virtual int GetNumPoints() const = 0;

    //! @brief ... return point id's
    //! @return ... array of point id's
    virtual const unsigned int* GetPoints() const = 0;

    //! @brief ... returns the corresponding Vtk cell type
    //! @return ... Vtk cell type
    virtual unsigned int GetVtkCellType() const = 0;

    //! @brief ... set tensor data
    //! @param data ... data
    void SetData(int dataIndex, Eigen::VectorXd data);

    //! @param dataIndex ... data index
    const Eigen::VectorXd& GetData(int dataIndex) const;

protected:

    std::vector<Eigen::VectorXd> mData;
};

}

