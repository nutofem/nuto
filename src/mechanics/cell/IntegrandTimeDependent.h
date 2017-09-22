#pragma once

namespace NuTo
{

#include "mechanics/nodes/DofMatrixContainer.h"

template <int TDim>
class CellData;

template <int TDim>
class CellIPData;


template <class T>
class DofVector;

template <int TDim>
class IntegrandTimeDependent
{
public:
    static constexpr int Dim = TDim;
    virtual NuTo::IntegrandTimeDependent<TDim>* Clone() const = 0;
    virtual NuTo::DofVector<double> Gradient(const NuTo::CellData<TDim>& rCellData, const NuTo::CellIPData<TDim>& rCellIPData) = 0;
    virtual NuTo::DofMatrix<double> Hessian0(const NuTo::CellData<TDim>& rCellData,
                                             const NuTo::CellIPData<TDim>& rCellIPData) const = 0;

    //virtual NuTo::DofMatrix<double> Hessian1(double t, double delta_t, const NuTo::CellData& rCellData,
    //                                         const NuTo::CellIPData& rCellIPData) const = 0;

    //virtual NuTo::DofMatrix<double> Hessian2(double t, double delta_t, const NuTo::CellData& rCellData,
    //                                         const NuTo::CellIPData& rCellIPData) const = 0;

    //virtual void UpdateHistoryData(double t, double delta_t, const NuTo::CellData& rCellData,
    //                               const NuTo::CellIPData& rCellIPData) = 0;

private:
};

} /* NuTo */
