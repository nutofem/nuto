#pragma once

#include <vector>
#include <eigen3/Eigen/Core>
#include "visualize/VisualizeEnum.h"

namespace NuTo
{

//! @brief ... base class for visualization cells
class CellBase
{
public:
    //! @brief constructor
    //! @param pointIds ... point ids
    //! @param cellType ... cell type enum
    //! @param numData ... number of different data fields
    CellBase(std::vector<int> pointIds, eCellTypes cellType, int numData);

    //! @brief ... return number of cell points
    //! @return ... number of cell points
    int GetNumPoints() const;

    //! @brief ... return point id's
    //! @return ... array of point id's
    const std::vector<int>& GetPointIds() const;

    //! @brief ... setter for mPointIds
    void SetPointIds(std::vector<int> pointIds);

    //! @brief ... returns the corresponding cell type
    //! @return ... cell type
    eCellTypes GetCellType() const;

    //! @brief ... set tensor data
    //! @param data ... data
    void SetData(int dataIndex, Eigen::VectorXd data);

    //! @param dataIndex ... data index
    const Eigen::VectorXd& GetData(int dataIndex) const;

protected:
    std::vector<Eigen::VectorXd> mData;
    std::vector<int> mPointIds;
    eCellTypes mCellType;
};
} // namespace NuTo
