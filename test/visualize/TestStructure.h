#pragma once

#include <fakeit.hpp>
#include "base/Group.h"
#include "mechanics/cell/CellInterface.h"
#include "mechanics/elements/ElementCollection.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"

namespace NuTo
{
namespace Test
{

class VisualizeTestStructure
{
    /* (0,4)      (6,5)   (8,5)
     *   3---------2-------5
     *   |         |       |
     *   |         |       |
     *   |         |       |
     *   |         |       |
     *   0---------1-------4
     * (1,1)      (5,0)   (8,1)
     */

    NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d(1, 1));
    NuTo::NodeSimple n1 = NuTo::NodeSimple(Eigen::Vector2d(5, 0));
    NuTo::NodeSimple n2 = NuTo::NodeSimple(Eigen::Vector2d(6, 5));
    NuTo::NodeSimple n3 = NuTo::NodeSimple(Eigen::Vector2d(0, 4));
    NuTo::NodeSimple n4 = NuTo::NodeSimple(Eigen::Vector2d(8, 1));
    NuTo::NodeSimple n5 = NuTo::NodeSimple(Eigen::Vector2d(8, 5));

    NuTo::NodeSimple nd0 = NuTo::NodeSimple(Eigen::Vector2d(1, 1) / 10);
    NuTo::NodeSimple nd1 = NuTo::NodeSimple(Eigen::Vector2d(5, 0) / 10);
    NuTo::NodeSimple nd2 = NuTo::NodeSimple(Eigen::Vector2d(6, 5) / 10);
    NuTo::NodeSimple nd3 = NuTo::NodeSimple(Eigen::Vector2d(0, 4) / 10);
    NuTo::NodeSimple nd4 = NuTo::NodeSimple(Eigen::Vector2d(8, 1) / 10);
    NuTo::NodeSimple nd5 = NuTo::NodeSimple(Eigen::Vector2d(8, 5) / 10);

    NuTo::InterpolationQuadLinear interpolation = NuTo::InterpolationQuadLinear();

    NuTo::ElementCollectionFem elements0 =
            NuTo::ElementCollectionFem(NuTo::ElementFem({n0, n1, n2, n3}, interpolation));
    NuTo::ElementCollectionFem elements1 =
            NuTo::ElementCollectionFem(NuTo::ElementFem({n1, n4, n5, n2}, interpolation));

    fakeit::Mock<NuTo::CellInterface> cell0;
    fakeit::Mock<NuTo::CellInterface> cell1;

public:
    VisualizeTestStructure(NuTo::DofType d)
    {
        elements0.AddDofElement(d, NuTo::ElementFem({nd0, nd1, nd2, nd3}, interpolation));
        elements1.AddDofElement(d, NuTo::ElementFem({nd1, nd4, nd5, nd2}, interpolation));
    }

    NuTo::Group<NuTo::CellInterface> Cells()
    {
        Eigen::VectorXd bla = Eigen::Vector3d::Constant(1.0);
        fakeit::When(ConstOverloadedMethod(cell0, Interpolate, Eigen::VectorXd(Eigen::VectorXd)))
                .AlwaysReturn(bla);
        fakeit::When(ConstOverloadedMethod(cell1, Interpolate, Eigen::VectorXd(Eigen::VectorXd)))
                .AlwaysReturn(bla);
        return {cell0.get(), cell1.get()};
    }
};

} /* Test */

} /* NuTo */
