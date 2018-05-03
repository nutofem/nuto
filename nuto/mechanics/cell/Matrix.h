#pragma once
#include "nuto/mechanics/interpolation/TypeDefs.h"

namespace NuTo
{
namespace Matrix
{
inline Eigen::MatrixXd N(const Eigen::VectorXd& shapeFunctions, int numNodes, int dim)
{
    Eigen::MatrixXd n(dim, dim * numNodes);
    for (int i = 0; i < numNodes; ++i)
        n.block(0, i * dim, dim, dim) = Eigen::MatrixXd::Identity(dim, dim) * shapeFunctions[i];
    return n;
}
} /*Matrix */
} /* NuTo */
