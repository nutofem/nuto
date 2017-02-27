#include "mechanics/interpolation/Interpolation.h"
#include "mechanics/MechanicsException.h"

Eigen::MatrixXd NuTo::Interpolation::GetN(const Eigen::VectorXd& rIPCoords) const
{
    int dim = GetDimension();
    Eigen::MatrixXd N(dim, dim * GetNumNodes());

    auto shapeFunctions = GetShapeFunctions(rIPCoords);

    for (size_t i = 0; i < GetNumNodes(); ++i)
        N.block(0, i * dim, dim, dim) = Eigen::MatrixXd::Identity(dim, dim) * shapeFunctions[i];
    return N;
}
Eigen::MatrixXd NuTo::Interpolation::GetBGradient(const Eigen::VectorXd& rIPCoords) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented.");
}
Eigen::MatrixXd NuTo::Interpolation::GetBStrain(const Eigen::VectorXd& rIPCoords) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented.");
}
