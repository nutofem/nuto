#pragma once

namespace NuTo
{

#include "mechanics/nodes/DofMatrixContainer.h"

class CellData;
class CellIPData;


template <class T>
class DofVector;

class IntegrandTimeDependent
{
public:

    virtual NuTo::IntegrandTimeDependent* Clone() const = 0;

    virtual NuTo::DofVector<double> Gradient(const NuTo::CellData& rCellData, const NuTo::CellIPData& rCellIPData) = 0;

    virtual NuTo::DofMatrix<double> Hessian0(const NuTo::CellData& rCellData,
                                             const NuTo::CellIPData& rCellIPData) const = 0;

    //virtual NuTo::DofMatrix<double> Hessian1(double t, double delta_t, const NuTo::CellData& rCellData,
    //                                         const NuTo::CellIPData& rCellIPData) const = 0;

    //virtual NuTo::DofMatrix<double> Hessian2(double t, double delta_t, const NuTo::CellData& rCellData,
    //                                         const NuTo::CellIPData& rCellIPData) const = 0;

    //virtual void UpdateHistoryData(double t, double delta_t, const NuTo::CellData& rCellData,
    //                               const NuTo::CellIPData& rCellIPData) = 0;

private:
};

} /* NuTo */
