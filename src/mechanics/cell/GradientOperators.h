#pragma once

#include <Eigen/Core>
#include "base/Exception.h"

namespace NuTo
{

namespace B
{
struct Interface
{
    virtual Eigen::MatrixXd operator()(const Eigen::MatrixXd& dNdX) const = 0;
};

struct Gradient : Interface
{
    Eigen::MatrixXd operator()(const Eigen::MatrixXd& dNdX) const override
    {
        return dNdX.transpose();
    }
};

struct Strain : Interface
{
    Eigen::MatrixXd operator()(const Eigen::MatrixXd& dNdX) const override
    {
        const int dim = dNdX.cols();
        const int numNodes = dNdX.rows();
        switch (dim)
        {
        case 1:
            return dNdX.transpose();
        case 2:
        {
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, numNodes * 2);
            for (int iNode = 0, iColumn = 0; iNode < numNodes; ++iNode, iColumn += 2)
            {
                double dNdXx = dNdX(iNode, 0);
                double dNdXy = dNdX(iNode, 1);

                B(0, iColumn) = dNdXx;
                B(1, iColumn + 1) = dNdXy;
                B(2, iColumn) = dNdXy;
                B(2, iColumn + 1) = dNdXx;
            }
            return B;
        }
        case 3:
        {
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, numNodes * 3);
            for (int iNode = 0, iColumn = 0; iNode < numNodes; ++iNode, iColumn += 3)
            {
                double dNdXx = dNdX(iNode, 0);
                double dNdXy = dNdX(iNode, 1);
                double dNdXz = dNdX(iNode, 2);

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


                B(0, iColumn) = dNdXx;
                B(1, iColumn + 1) = dNdXy;
                B(2, iColumn + 2) = dNdXz;

                B(3, iColumn + 1) = dNdXz;
                B(3, iColumn + 2) = dNdXy;

                B(4, iColumn) = dNdXz;
                B(4, iColumn + 2) = dNdXx;

                B(5, iColumn) = dNdXy;
                B(5, iColumn + 1) = dNdXx;
            }
            return B;
        }
        default:
            throw Exception(__PRETTY_FUNCTION__, "c'mon.");
        }
    }
};
} /* B */
} /* NuTo */
