#pragma once
#include "mechanics/interpolation/TypeDefs.h"

namespace NuTo
{
namespace Matrix
{
inline NMatrix N(const ShapeFunctions& shapeFunctions, int numNodes, int dim)
{
    NMatrix n(dim, dim * numNodes);
    for (int i = 0; i < numNodes; ++i)
        n.block(0, i * dim, dim, dim) = Eigen::MatrixXd::Identity(dim, dim) * shapeFunctions[i];
    return n;
}
} /*Matrix */
} /* NuTo */
