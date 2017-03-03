#include "mechanics/interpolation/Interpolation.h"
#include "mechanics/MechanicsException.h"

NuTo::NMatrix NuTo::Interpolation::GetN(const Eigen::VectorXd& rLocalIPCoords) const
{
    int dim = GetDofDimension();
    Eigen::MatrixXd N(dim, dim * GetNumNodes());

    auto shapeFunctions = GetShapeFunctions(rLocalIPCoords);

    for (size_t i = 0; i < GetNumNodes(); ++i)
        N.block(0, i * dim, dim, dim) = Eigen::MatrixXd::Identity(dim, dim) * shapeFunctions[i];
    return N;
}

NuTo::BGradientLocal NuTo::Interpolation::GetBGradient(const Eigen::VectorXd& rLocalIPCoords) const
{
    if (GetDofDimension() != 1)
        throw MechanicsException(__PRETTY_FUNCTION__, "moep. Only for scalars.");
    return GetDerivativeShapeFunctions(rLocalIPCoords).transpose();
}

NuTo::BStrainLocal NuTo::Interpolation::GetBStrain(const Eigen::VectorXd& rLocalIPCoords) const
{
    auto derivativeShapeFunctions = GetDerivativeShapeFunctions(rLocalIPCoords);
    switch (derivativeShapeFunctions.cols())
    {
    case 1:
        return derivativeShapeFunctions.transpose();
    case 2:
    {
        BStrainLocal B = Eigen::MatrixXd::Zero(3, GetNumNodes() * 2);
        for (int iNode = 0, iColumn = 0; iNode < GetNumNodes(); ++iNode, iColumn += 2)
        {
            double dNdX = derivativeShapeFunctions(iNode, 0);
            double dNdY = derivativeShapeFunctions(iNode, 1);

            B(0, iColumn)     = dNdX;
            B(1, iColumn + 1) = dNdY;
            B(2, iColumn)     = dNdY;
            B(2, iColumn + 1) = dNdX;
        }
        return B;
    }
    case 3:
    {

        BStrainLocal B = Eigen::MatrixXd::Zero(6, GetNumNodes() * 3);

        for (int iNode = 0, iColumn = 0; iNode < GetNumNodes(); ++iNode, iColumn += 3)
        {
            double dNdX = derivativeShapeFunctions(iNode, 0);
            double dNdY = derivativeShapeFunctions(iNode, 1);
            double dNdZ = derivativeShapeFunctions(iNode, 2);

            /* according to Jirásek
             *
             *     +0  +1  +2
             *    -------------
             * 0 |  dx  0   0  |     - e_x
             * 1 |  0   dy  0  |     - e_y
             * 2 |  0   0   dz |     - e_z
             * 3 |  0   dz  dy |     - g_yz
             * 4 |  dz  0   dx |     - g_xz
             * 5 |  dy  dx  0  |     - g_xy
             *    -------------
             */


            B(0, iColumn)     = dNdX;
            B(1, iColumn + 1) = dNdY;
            B(2, iColumn + 2) = dNdZ;

            B(3, iColumn + 1) = dNdZ;
            B(3, iColumn + 2) = dNdY;

            B(4, iColumn)     = dNdZ;
            B(4, iColumn + 2) = dNdX;

            B(5, iColumn)     = dNdY;
            B(5, iColumn + 1) = dNdX;
        }
        return B;
    }
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "c'mon.");
    }
}
