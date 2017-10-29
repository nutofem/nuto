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

    NuTo::InterpolationQuadLinear interpolation = NuTo::InterpolationQuadLinear(2);

    NuTo::ElementCollectionFem elements0 =
            NuTo::ElementCollectionFem(NuTo::ElementFem({n0, n1, n2, n3}, interpolation));
    NuTo::ElementCollectionFem elements1 =
            NuTo::ElementCollectionFem(NuTo::ElementFem({n1, n4, n5, n2}, interpolation));

    NuTo::IpValues ip00 = {{Eigen::Vector3d(1, 1, 1), "Stress"}, {Eigen::Vector3d(10, 10, 10), "Strain"}};
    NuTo::IpValues ip01 = {{Eigen::Vector3d(2, 2, 2), "Stress"}, {Eigen::Vector3d(20, 20, 20), "Strain"}};

    NuTo::IpValues ip10 = {{Eigen::Vector3d(3, 3, 3), "Stress"}, {Eigen::Vector3d(30, 30, 30), "Strain"}};
    NuTo::IpValues ip11 = {{Eigen::Vector3d(4, 4, 4), "Stress"}, {Eigen::Vector3d(40, 40, 40), "Strain"}};

    fakeit::Mock<NuTo::CellInterface> cell0;
    fakeit::Mock<NuTo::CellInterface> cell1;

public:
    VisualizeTestStructure(NuTo::DofType d)
    {
        elements0.AddDofElement(d, NuTo::ElementFem({nd0, nd1, nd2, nd3}, interpolation));
        elements1.AddDofElement(d, NuTo::ElementFem({nd1, nd4, nd5, nd2}, interpolation));
    }

    NuTo::Groups::Group<NuTo::CellInterface> Cells()
    {

        Method(cell0, GetElementCollection) = elements0;
        Method(cell0, GetIpValues) = {ip00, ip01};
        Method(cell1, GetElementCollection) = elements1;
        Method(cell1, GetIpValues) = {ip10, ip11};

        return {cell0.get(), cell1.get()};
    }
};

} /* Test */

} /* NuTo */
