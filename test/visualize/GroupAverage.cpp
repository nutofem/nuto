#include "BoostUnitTest.h"
#include <fakeit.hpp>

#include "mechanics/elements/ElementCollection.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "visualize/GroupAverage.h"


BOOST_AUTO_TEST_CASE(Visualize)
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
    NuTo::NodeSimple n0(Eigen::Vector2d(1, 1));
    NuTo::NodeSimple n1(Eigen::Vector2d(5, 0));
    NuTo::NodeSimple n2(Eigen::Vector2d(6, 5));
    NuTo::NodeSimple n3(Eigen::Vector2d(0, 4));
    NuTo::NodeSimple n4(Eigen::Vector2d(8, 1));
    NuTo::NodeSimple n5(Eigen::Vector2d(8, 5));

    NuTo::InterpolationQuadLinear interpolation(2);

    NuTo::ElementCollectionFem elements0(NuTo::ElementFem({n0, n1, n2, n3}, interpolation));
    NuTo::ElementCollectionFem elements1(NuTo::ElementFem({n1, n4, n5, n2}, interpolation));

    NuTo::DofType dof("NodeCoordinatesDiv10", 2);

    NuTo::NodeSimple nd0(Eigen::Vector2d(1, 1) / 10);
    NuTo::NodeSimple nd1(Eigen::Vector2d(5, 0) / 10);
    NuTo::NodeSimple nd2(Eigen::Vector2d(6, 5) / 10);
    NuTo::NodeSimple nd3(Eigen::Vector2d(0, 4) / 10);
    NuTo::NodeSimple nd4(Eigen::Vector2d(8, 1) / 10);
    NuTo::NodeSimple nd5(Eigen::Vector2d(8, 5) / 10);

    elements0.AddDofElement(dof, NuTo::ElementFem({nd0, nd1, nd2, nd3}, interpolation));
    elements1.AddDofElement(dof, NuTo::ElementFem({nd1, nd4, nd5, nd2}, interpolation));

    fakeit::Mock<NuTo::CellInterface> cell0;
    fakeit::Mock<NuTo::CellInterface> cell1;
    Method(cell0, GetElementCollection) = elements0;
    Method(cell1, GetElementCollection) = elements1;


    /*
     *  Actual visualize stuff begins here.
     */

    NuTo::Visualize::CellGeometry cellGeometry;
    cellGeometry.mCornerCoords = {interpolation.GetLocalCoords(0), interpolation.GetLocalCoords(1),
                                  interpolation.GetLocalCoords(2), interpolation.GetLocalCoords(3)};
    cellGeometry.mCellType = NuTo::eCellTypes::QUAD;


    NuTo::Visualize::GroupAverage v({cell0.get(), cell1.get()}, cellGeometry);
    v.ExtractGeometry();
    v.Visualize("PdeAverage.vtu", {dof}, false);
}
