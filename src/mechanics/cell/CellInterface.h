#pragma once

#include "mechanics/cell/IpValues.h"
#include "mechanics/dofs/DofVector.h"
#include "mechanics/dofs/DofMatrix.h"
#include <functional>

namespace NuTo
{
class CellData;
class CellIpData;

class CellInterface
{
public:
    virtual ~CellInterface() = default;

    using ScalarFunction = std::function<double(const CellData&, const CellIpData&)>;
    using VectorFunction = std::function<DofVector<double>(const CellData&, const CellIpData&)>;
    using MatrixFunction = std::function<DofMatrix<double>(const CellData&, const CellIpData&)>;

    using VoidFunction = std::function<void(const CellData&, const CellIpData&)>;

    virtual double Integrate(ScalarFunction) = 0;
    virtual DofVector<double> Integrate(VectorFunction) = 0;
    virtual DofMatrix<double> Integrate(MatrixFunction) = 0;
    virtual void Apply(VoidFunction) = 0;

    virtual Eigen::VectorXi DofNumbering(DofType dof) = 0;
};
} /* NuTo */
