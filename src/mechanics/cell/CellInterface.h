#pragma once

#include "mechanics/cell/IpValues.h"
#include "mechanics/dofs/DofVector.h"
#include "mechanics/dofs/DofMatrix.h"
#include "mechanics/integrands/Operations.h"

namespace NuTo
{

class ElementCollection;

class CellInterface
{
public:
    virtual ~CellInterface() = default;

    virtual DofVector<double> Integrate(const VectorOperation&) = 0;
    virtual DofMatrix<double> Integrate(const MatrixOperation&) = 0;
    virtual double Integrate(const ScalarOperation&) = 0;
    virtual void Apply(const VoidOperation&) = 0;
    virtual Eigen::VectorXi DofNumbering(DofType dof) = 0;

    virtual std::vector<IpValues> GetIpValues() = 0;

    virtual std::vector<Eigen::VectorXd> GetIpCoordinates() = 0;

    virtual const ElementCollection& GetElementCollection() const = 0;
};
} /* NuTo */
