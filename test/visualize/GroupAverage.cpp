#include "BoostUnitTest.h"
#include <fakeit.hpp>

#include "mechanics/elements/ElementCollection.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "visualize/GroupAverage.h"


BOOST_AUTO_TEST_CASE(Visualize)
{
    /* (0,4)      (6,5)
     *   3---------2
     *   |         |
     *   |         |
     *   |         |
     *   |         |
     *   0---------1
     * (1,1)      (5,0)
     */
    NuTo::NodeSimple n0(Eigen::Vector2d(1, 1));
    NuTo::NodeSimple n1(Eigen::Vector2d(5, 0));
    NuTo::NodeSimple n2(Eigen::Vector2d(6, 5));
    NuTo::NodeSimple n3(Eigen::Vector2d(0, 4));

    NuTo::InterpolationQuadLinear interpolation(2);

    NuTo::ElementCollectionFem elements(NuTo::ElementFem({n0, n1, n2, n3}, interpolation));

    NuTo::DofType dof("NodeCoordinatesDiv10", 2);

    NuTo::NodeSimple nd0(Eigen::Vector2d(1, 1) / 10);
    NuTo::NodeSimple nd1(Eigen::Vector2d(5, 0) / 10);
    NuTo::NodeSimple nd2(Eigen::Vector2d(6, 5) / 10);
    NuTo::NodeSimple nd3(Eigen::Vector2d(0, 4) / 10);

    elements.AddDofElement(dof, NuTo::ElementFem({nd0, nd1, nd2, nd3}, interpolation));

    fakeit::Mock<NuTo::CellInterface> cell;
    Method(cell, GetElementCollection) = elements;

    /*
     *  Actual visualize stuff begins here.
     */

    NuTo::Visualize::CellGeometry cellGeometry;
    cellGeometry.mCornerCoords = {interpolation.GetLocalCoords(0), interpolation.GetLocalCoords(1),
                                  interpolation.GetLocalCoords(2), interpolation.GetLocalCoords(3)};
    cellGeometry.mCellType = NuTo::eCellTypes::QUAD;


    NuTo::Visualize::GroupAverage v({cell.get()}, cellGeometry);
    v.ExtractGeometry();
    v.Visualize("PdeAverage.vtu", {dof}, false);
}
