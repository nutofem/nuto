#pragma once

#include <vector>
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

    using ScalarFunction = std::function<double(const CellIpData&)>;
    using VectorFunction = std::function<DofVector<double>(const CellIpData&)>;
    using MatrixFunction = std::function<DofMatrix<double>(const CellIpData&)>;

    using VoidFunction = std::function<void(const CellIpData&)>;
    using EvalFunction = std::function<Eigen::VectorXd(const CellIpData&)>;

    virtual double Integrate(ScalarFunction) = 0;
    virtual DofVector<double> Integrate(VectorFunction) = 0;
    virtual DofMatrix<double> Integrate(MatrixFunction) = 0;
    virtual void Apply(VoidFunction) = 0;

    virtual std::vector<Eigen::VectorXd> Eval(EvalFunction f) const = 0;

    virtual Eigen::VectorXi DofNumbering(DofType dof) = 0;

    //! Coordinate interpolation
    virtual Eigen::VectorXd Interpolate(Eigen::VectorXd naturalCoords) const = 0;
    //! Dof interpolation
    virtual Eigen::VectorXd Interpolate(Eigen::VectorXd naturalCoords, DofType dof) const = 0;
};
} /* NuTo */
