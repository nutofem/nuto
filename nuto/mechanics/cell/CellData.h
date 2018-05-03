#pragma once

#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/elements/ElementCollection.h"

namespace NuTo
{


//! @brief Extracts 'cell data' like nodal values from the cell
//! @remark The life time of objects of this class is limited to the evaluation of one NuTo::Cell. This class will
//! calculate (and cache) all the information that are required once per cell. Classic example: NodeValues. Even if
//! `GetNodeValues()` is called 42 times, they should only be extracted once.
class CellData
{
public:
    CellData(const ElementCollection& elements, int cellId)
        : mElements(elements)
        , mCellId(cellId)
    {
    }

    int GetCellId() const
    {
        return mCellId;
    }

    Eigen::VectorXd GetCoordinates() const
    {
        return mElements.CoordinateElement().ExtractNodeValues();
    }

    const Eigen::VectorXd& GetNodeValues(DofType dofType, int instance = 0) const
    {
        if (instance >= static_cast<int>(mNodeValues[dofType].size()))
            mNodeValues[dofType].resize(instance + 1);

        Eigen::VectorXd& nodeValues = mNodeValues[dofType][instance];
        if (nodeValues.size() == 0)
            nodeValues = mElements.DofElement(dofType).ExtractNodeValues(instance);

        return nodeValues;
    }

    const ElementCollection& Elements() const
    {
        return mElements;
    }

private:
    mutable DofContainer<std::vector<Eigen::VectorXd>> mNodeValues;
    const ElementCollection& mElements;
    int mCellId;
};
} /* NuTo */
