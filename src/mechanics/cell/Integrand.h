#pragma once
#include <memory>
#include "mechanics/cell/CellData.h"
#include "mechanics/cell/CellIPData.h"
#include "mechanics/nodes/DofVector.h"

namespace NuTo
{
struct IPValue
{
    std::string mName;
    Eigen::MatrixXd mValue;
};

class Integrand
{
public:
    virtual std::unique_ptr<Integrand> Clone() const = 0;
    virtual ~Integrand()        = default;
    virtual DofVector Gradient(const CellData&, const CellIPData&) = 0;
    virtual std::vector<IPValue> IPValues(const CellData&, const CellIPData&) = 0;
};
} /* NuTo */
