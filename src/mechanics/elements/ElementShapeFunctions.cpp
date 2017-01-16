/*
 * ElementShapeFunctions.cpp
 *
 *  Created on: 30 Mar 2015
 *      Author: ttitsche
 */
#include <assert.h>
#include "mechanics/elements/ElementShapeFunctions.h"
#include "mechanics/MechanicsException.h"

#include <iostream>

namespace NuTo
{
namespace ShapeFunctions1D // interval -1 to 1
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1: return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..1)");
    }
}

Eigen::Matrix<double, 2, 1> ShapeFunctionsTrussOrder1(const Eigen::VectorXd& rCoordinates)
{
    return Eigen::Vector2d(0.5 * (1. - rCoordinates[0]),
                           0.5 * (1. + rCoordinates[0]));
}

Eigen::Matrix<double, 2, 1> DerivativeShapeFunctionsTrussOrder1(const Eigen::VectorXd& rCoordinates)
{
    return Eigen::Vector2d(-0.5, 0.5);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1: return Eigen::Matrix<double, 1, 1>::Constant(0.);
    case 2: return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..2)");
    }
}

Eigen::Matrix<double, 3, 1> ShapeFunctionsTrussOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 3, 1> shapeFunctions;
    shapeFunctions[0] = 0.5 * (1. - rCoordinates[0]) - 0.5 * (1. - rCoordinates[0] * rCoordinates[0]);
    shapeFunctions[1] = 1. - rCoordinates[0] * rCoordinates[0];
    shapeFunctions[2] = 0.5 * (1. + rCoordinates[0]) - 0.5 * (1. - rCoordinates[0] * rCoordinates[0]);
    return shapeFunctions;
}

Eigen::Matrix<double, 3, 1> DerivativeShapeFunctionsTrussOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 3, 1> derivativeShapeFunctions;
    derivativeShapeFunctions[0] = -0.5 + rCoordinates[0];
    derivativeShapeFunctions[1] = -2.0 * rCoordinates[0];
    derivativeShapeFunctions[2] = 0.5 + rCoordinates[0];
    return derivativeShapeFunctions;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder3(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1: return Eigen::Matrix<double, 1, 1>::Constant(-1. / 3.);
    case 2: return Eigen::Matrix<double, 1, 1>::Constant(1. / 3.);
    case 3: return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..3)");
    }
}

Eigen::Matrix<double, 4, 1> ShapeFunctionsTrussOrder3(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> shapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    double r3 = r2 * r;
    shapeFunctions[0] = -0.0625 + 0.0625 * r + 0.5625 * r2 - 0.5625 * r3;
    shapeFunctions[1] = +0.5625 - 1.6875 * r - 0.5625 * r2 + 1.6875 * r3;
    shapeFunctions[2] = +0.5625 + 1.6875 * r - 0.5625 * r2 - 1.6875 * r3;
    shapeFunctions[3] = -0.0625 - 0.0625 * r + 0.5625 * r2 + 0.5625 * r3;
    return shapeFunctions;

}

Eigen::Matrix<double, 4, 1> DerivativeShapeFunctionsTrussOrder3(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    derivativeShapeFunctions[0] = +0.0625 + 1.125 * r - 1.6875 * r2;
    derivativeShapeFunctions[1] = -1.6875 - 1.125 * r + 5.0625 * r2;
    derivativeShapeFunctions[2] = +1.6875 - 1.125 * r - 5.0625 * r2;
    derivativeShapeFunctions[3] = -0.0625 + 1.125 * r + 1.6875 * r2;
    return derivativeShapeFunctions;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder4(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1: return Eigen::Matrix<double, 1, 1>::Constant(-0.5);
    case 2: return Eigen::Matrix<double, 1, 1>::Constant(0.);
    case 3: return Eigen::Matrix<double, 1, 1>::Constant(0.5);
    case 4: return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..4)");
    }
}

Eigen::Matrix<double, 5, 1> ShapeFunctionsTrussOrder4(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 5, 1> shapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    double r3 = r2 * r;
    double r4 = r3 * r;
    shapeFunctions[0] = +0.166666666667 * r - 0.166666666667 * r2 - 0.666666666667 * r3 + 0.666666666667 * r4;
    shapeFunctions[1] = -1.33333333333 * r + 2.66666666667 * r2 + 1.33333333333 * r3 - 2.66666666667 * r4;
    shapeFunctions[2] = +1.0 - 5.0 * r2 + 4.0 * r4;
    shapeFunctions[3] = +1.33333333333 * r + 2.66666666667 * r2 - 1.33333333333 * r3 - 2.66666666667 * r4;
    shapeFunctions[4] = -0.166666666667 * r - 0.166666666667 * r2 + 0.666666666667 * r3 + 0.666666666667 * r4;
    return shapeFunctions;
}

Eigen::Matrix<double, 5, 1> DerivativeShapeFunctionsTrussOrder4(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 5, 1> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    double r3 = r2 * r;
    derivativeShapeFunctions[0] = +0.166666666667 - 0.333333333333 * r - 2.0 * r2 + 2.66666666667 * r3;
    derivativeShapeFunctions[1] = -1.33333333333 + 5.33333333333 * r + 4.0 * r2 - 10.6666666667 * r3;
    derivativeShapeFunctions[2] = -10.0 * r + 16.0 * r3;
    derivativeShapeFunctions[3] = +1.33333333333 + 5.33333333333 * r - 4.0 * r2 - 10.6666666667 * r3;
    derivativeShapeFunctions[4] = -0.166666666667 - 0.333333333333 * r + 2.0 * r2 + 2.66666666667 * r3;
    return derivativeShapeFunctions;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussSpectralOrder3(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1: return Eigen::Matrix<double, 1, 1>::Constant(-0.447213595499957928);
    case 2: return Eigen::Matrix<double, 1, 1>::Constant(0.447213595499957928);
    case 3: return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..3)");
    }
}

Eigen::Matrix<double, 4, 1> ShapeFunctionsTrussSpectralOrder3(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> shapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    double r3 = r * r2;
    shapeFunctions[0] = -0.125 + 0.125 * r + 0.625 * r2 - 0.625 * r3;
    shapeFunctions[1] = 0.625 - 1.3975424859373684 * r - 0.625 * r2 + 1.3975424859373684 * r3;
    shapeFunctions[2] = 0.625 + 1.3975424859373684 * r - 0.625 * r2 - 1.3975424859373684 * r3;
    shapeFunctions[3] = -0.125 - 0.125 * r + 0.625 * r2 + 0.625 * r3;
    return shapeFunctions;
}

Eigen::Matrix<double, 4, 1> DerivativeShapeFunctionsTrussSpectralOrder3(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    derivativeShapeFunctions[0] = 0.125 + 1.25 * r - 1.875 * r2;
    derivativeShapeFunctions[1] = -1.3975424859373684 - 1.25 * r + 4.192627457812105 * r2;
    derivativeShapeFunctions[2] = 1.3975424859373684 - 1.25 * r - 4.192627457812105 * r2;
    derivativeShapeFunctions[3] = -0.125 + 1.25 * r + 1.875 * r2;
    return derivativeShapeFunctions;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussSpectralOrder4(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Matrix<double, 1, 1>::Constant(-1.);
    case 1: return Eigen::Matrix<double, 1, 1>::Constant(-0.654653670707977087);
    case 2: return Eigen::Matrix<double, 1, 1>::Constant(0.0);
    case 3: return Eigen::Matrix<double, 1, 1>::Constant(0.654653670707977087);
    case 4: return Eigen::Matrix<double, 1, 1>::Constant(1.);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..4)");
    }
}

Eigen::Matrix<double, 5, 1> ShapeFunctionsTrussSpectralOrder4(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 5, 1> shapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    double r3 = r2 * r;
    double r4 = r3 * r;
    shapeFunctions[0] = +0.375 * r - 0.375 * r2 - 0.875 * r3 + 0.875 * r4;
    shapeFunctions[1] = -1.336584577695453 * r + 2.041666666666666 * r2 + 1.336584577695453 * r3 - 2.041666666666666 * r4;
    shapeFunctions[2] = 1. - 3.333333333333333 * r2 + 2.333333333333333 * r4;
    shapeFunctions[3] = +1.336584577695453 * r + 2.041666666666666 * r2 - 1.336584577695453 * r3 - 2.041666666666666 * r4;
    shapeFunctions[4] = -0.375 * r - 0.375 * r2 + 0.875 * r3 + 0.875 * r4;
    return shapeFunctions;
}

Eigen::Matrix<double, 5, 1> DerivativeShapeFunctionsTrussSpectralOrder4(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 5, 1> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double r2 = r * r;
    double r3 = r2 * r;
    derivativeShapeFunctions[0] = 0.375 - 0.75 * r - 2.625 * r2 + 3.5 * r3;
    derivativeShapeFunctions[1] = -1.336584577695453 + 4.083333333333333 * r + 4.009753733086359517 * r2 - 8.16666666666666 * r3;
    derivativeShapeFunctions[2] = -6.666666666666666 * r + 9.33333333333333 * r3;
    derivativeShapeFunctions[3] = 1.336584577695453 + 4.083333333333333 * r - 4.009753733086359517 * r2 - 8.16666666666666 * r3;
    derivativeShapeFunctions[4] = -0.375 - 0.75 * r + 2.625 * r2 + 3.5 * r3;
    return derivativeShapeFunctions;
}

}

namespace ShapeFunctions2D
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Vector2d(0.0, 0.0);
    case 1: return Eigen::Vector2d(1.0, 0.0);
    case 2: return Eigen::Vector2d(0.0, 1.0);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..2)");
    }
}

Eigen::Matrix<double, 3, 1> ShapeFunctionsTriangleOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 3, 1> shapeFunctions;
    shapeFunctions[0] = 1. - rCoordinates(0) - rCoordinates(1);
    shapeFunctions[1] = rCoordinates(0);
    shapeFunctions[2] = rCoordinates(1);
    return shapeFunctions;
}

Eigen::Matrix<double, 3, 2> DerivativeShapeFunctionsTriangleOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 3, 2> derivativeShapeFunctions;
    derivativeShapeFunctions(0, 0) = -1.0;
    derivativeShapeFunctions(0, 1) = -1.0;

    derivativeShapeFunctions(1, 0) = 1.0;
    derivativeShapeFunctions(1, 1) = 0.0;

    derivativeShapeFunctions(2, 0) = 0.0;
    derivativeShapeFunctions(2, 1) = 1.0;
    return derivativeShapeFunctions;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Vector2d(0.0, 0.0);
    case 1: return Eigen::Vector2d(1.0, 0.0);
    case 2: return Eigen::Vector2d(0.0, 1.0);
    case 3: return Eigen::Vector2d(0.5, 0.0);
    case 4: return Eigen::Vector2d(0.5, 0.5);
    case 5: return Eigen::Vector2d(0.0, 0.5);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..5)");
    }
}

Eigen::Matrix<double, 6, 1> ShapeFunctionsTriangleOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 6, 1> shapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];

    shapeFunctions[0] = 2. * (r * r + s * s) + 4. * r * s - 3. * (r + s) + 1.;
    shapeFunctions[1] = 2. * r * r - r;
    shapeFunctions[2] = 2. * s * s - s;
    shapeFunctions[3] = -4. * r * (r + s - 1.);
    shapeFunctions[4] = 4. * r * s;
    shapeFunctions[5] = -4. * s * (s + r - 1.);
    return shapeFunctions;
}

Eigen::Matrix<double, 6, 2> DerivativeShapeFunctionsTriangleOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 6, 2> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];

    derivativeShapeFunctions(0, 0) = 4. * (r + s) - 3.;
    derivativeShapeFunctions(0, 1) = 4. * (r + s) - 3.;

    derivativeShapeFunctions(1, 0) = 4. * r - 1.;
    derivativeShapeFunctions(1, 1) = 0.;

    derivativeShapeFunctions(2, 0) = 0.;
    derivativeShapeFunctions(2, 1) = 4. * s - 1.;

    derivativeShapeFunctions(3, 0) = -8. * r - 4. * s + 4.;
    derivativeShapeFunctions(3, 1) = -4. * r;

    derivativeShapeFunctions(4, 0) = 4. * s;
    derivativeShapeFunctions(4, 1) = 4. * r;

    derivativeShapeFunctions(5, 0) = -4. * s;
    derivativeShapeFunctions(5, 1) = -8. * s - 4. * r + 4.;
    return derivativeShapeFunctions;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder3(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Vector2d(0.0, 0.0);
    case 1: return Eigen::Vector2d(1. / 3., 0.0);
    case 2: return Eigen::Vector2d(2. / 3., 0.0);
    case 3: return Eigen::Vector2d(1.0, 0.0);
    case 4: return Eigen::Vector2d(0.0, 1. / 3.);
    case 5: return Eigen::Vector2d(1. / 3., 1. / 3.);
    case 6: return Eigen::Vector2d(2. / 3., 1. / 3.);
    case 7: return Eigen::Vector2d(0., 2. / 3.);
    case 8: return Eigen::Vector2d(1. / 3., 2. / 3.);
    case 9: return Eigen::Vector2d(0., 1.);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..9)");
    }
}

Eigen::Matrix<double, 10, 1> ShapeFunctionsTriangleOrder3(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 10, 1> shapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];

    shapeFunctions[0] = +1.0 - 5.5 * r - 5.5 * s + 9.0 * r * r + 18.0 * r * s + 9.0 * s * s - 4.5 * r * r * r - 13.5 * r * r * s - 13.5 * r * s * s - 4.5 * s * s * s;
    shapeFunctions[1] = +9.0 * r - 22.5 * r * r - 22.5 * r * s + 13.5 * r * r * r + 27.0 * r * r * s + 13.5 * r * s * s;
    shapeFunctions[2] = -4.5 * r + 18.0 * r * r + 4.5 * r * s - 13.5 * r * r * r - 13.5 * r * r * s;
    shapeFunctions[3] = +1.0 * r - 4.5 * r * r + 4.5 * r * r * r;
    shapeFunctions[4] = +9.0 * s - 22.5 * r * s - 22.5 * s * s + 13.5 * r * r * s + 27.0 * r * s * s + 13.5 * s * s * s;
    shapeFunctions[5] = +27.0 * r * s - 27.0 * r * r * s - 27.0 * r * s * s;
    shapeFunctions[6] = -4.5 * r * s + 13.5 * r * r * s;
    shapeFunctions[7] = -4.5 * s + 4.5 * r * s + 18.0 * s * s - 13.5 * r * s * s - 13.5 * s * s * s;
    shapeFunctions[8] = -4.5 * r * s + 13.5 * r * s * s;
    shapeFunctions[9] = +1.0 * s - 4.5 * s * s + 4.5 * s * s * s;

    return shapeFunctions;
}

Eigen::Matrix<double, 10, 2> DerivativeShapeFunctionsTriangleOrder3(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 10, 2> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];

    derivativeShapeFunctions(0, 0) = -5.5 + 18.0 * r + 18.0 * s - 13.5 * r * r - 27.0 * r * s - 13.5 * s * s;
    derivativeShapeFunctions(0, 1) = -5.5 + 18.0 * r + 18.0 * s - 13.5 * r * r - 27.0 * r * s - 13.5 * s * s;

    derivativeShapeFunctions(1, 0) = +9.0 - 45.0 * r - 22.5 * s + 40.5 * r * r + 54.0 * r * s + 13.5 * s * s;
    derivativeShapeFunctions(1, 1) = -22.5 * r + 27.0 * r * r + 27.0 * r * s;

    derivativeShapeFunctions(2, 0) = -4.5 + 36.0 * r + 4.5 * s - 40.5 * r * r - 27.0 * r * s;
    derivativeShapeFunctions(2, 1) = +4.5 * r - 13.5 * r * r;

    derivativeShapeFunctions(3, 0) = +1.0 - 9.0 * r + 13.5 * r * r;
    derivativeShapeFunctions(3, 1) = 0.;

    derivativeShapeFunctions(4, 0) = -22.5 * s + 27.0 * r * s + 27.0 * s * s;
    derivativeShapeFunctions(4, 1) = +9.0 - 22.5 * r - 45.0 * s + 13.5 * r * r + 54.0 * r * s + 40.5 * s * s;

    derivativeShapeFunctions(5, 0) = +27.0 * s - 54.0 * r * s - 27.0 * s * s;
    derivativeShapeFunctions(5, 1) = +27.0 * r - 27.0 * r * r - 54.0 * r * s;

    derivativeShapeFunctions(6, 0) = -4.5 * s + 27.0 * r * s;
    derivativeShapeFunctions(6, 1) = -4.5 * r + 13.5 * r * r;

    derivativeShapeFunctions(7, 0) = +4.5 * s - 13.5 * s * s;
    derivativeShapeFunctions(7, 1) = -4.5 + 4.5 * r + 36.0 * s - 27.0 * r * s - 40.5 * s * s;

    derivativeShapeFunctions(8, 0) = -4.5 * s + 13.5 * s * s;
    derivativeShapeFunctions(8, 1) = -4.5 * r + 27.0 * r * s;

    derivativeShapeFunctions(9, 0) = 0.;
    derivativeShapeFunctions(9, 1) = +1.0 - 9.0 * s + 13.5 * s * s;

    return derivativeShapeFunctions;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder4(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:  return Eigen::Vector2d(.00, .00);
    case 1:  return Eigen::Vector2d(.25, .00);
    case 2:  return Eigen::Vector2d(.50, .00);
    case 3:  return Eigen::Vector2d(.75, .00);
    case 4:  return Eigen::Vector2d(1.00, .00);
    case 5:  return Eigen::Vector2d(.00, .25);
    case 6:  return Eigen::Vector2d(.25, .25);
    case 7:  return Eigen::Vector2d(.50, .25);
    case 8:  return Eigen::Vector2d(.75, .25);
    case 9:  return Eigen::Vector2d(.00, .50);
    case 10: return Eigen::Vector2d(.25, .50);
    case 11: return Eigen::Vector2d(.50, .50);
    case 12: return Eigen::Vector2d(.00, .75);
    case 13: return Eigen::Vector2d(.25, .75);
    case 14: return Eigen::Vector2d(.00, 1.00);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..14)");
    }
}

Eigen::Matrix<double, 15, 1> ShapeFunctionsTriangleOrder4(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 15, 1> shapeFunctions;
    double r(rCoordinates(0));
    double s(rCoordinates(1));

    shapeFunctions[0] = +1.0 - 8.33333333333 * r - 8.33333333333 * s + 23.3333333333 * r * r + 46.6666666667 * r * s + 23.3333333333 * s * s - 26.6666666667 * r * r * r - 80.0 * r * r * s - 80.0 * r * s * s - 26.6666666667 * s * s * s + 10.6666666667 * s * s * s * s + 42.6666666667 * r * s * s * s
            + 64.0 * r * r * s * s + 42.6666666667 * r * r * r * s + 10.6666666667 * r * r * r * r;
    shapeFunctions[1] = +16.0 * r - 69.3333333333 * r * r - 69.3333333333 * r * s + 96.0 * r * r * r + 192.0 * r * r * s + 96.0 * r * s * s - 42.6666666667 * r * s * s * s - 128.0 * r * r * s * s - 128.0 * r * r * r * s - 42.6666666667 * r * r * r * r;
    shapeFunctions[2] = -12.0 * r + 76.0 * r * r + 28.0 * r * s - 128.0 * r * r * r - 144.0 * r * r * s - 16.0 * r * s * s + 64.0 * r * r * s * s + 128.0 * r * r * r * s + 64.0 * r * r * r * r;
    shapeFunctions[3] = +5.33333333333 * r - 37.3333333333 * r * r - 5.33333333333 * r * s + 74.6666666667 * r * r * r + 32.0 * r * r * s - 42.6666666667 * r * r * r * s - 42.6666666667 * r * r * r * r;
    shapeFunctions[4] = -1.0 * r + 7.33333333333 * r * r - 16.0 * r * r * r + 10.6666666667 * r * r * r * r;
    shapeFunctions[5] = +16.0 * s - 69.3333333333 * r * s - 69.3333333333 * s * s + 96.0 * r * r * s + 192.0 * r * s * s + 96.0 * s * s * s - 42.6666666667 * s * s * s * s - 128.0 * r * s * s * s - 128.0 * r * r * s * s - 42.6666666667 * r * r * r * s;
    shapeFunctions[6] = +96.0 * r * s - 224.0 * r * r * s - 224.0 * r * s * s + 128.0 * r * s * s * s + 256.0 * r * r * s * s + 128.0 * r * r * r * s;
    shapeFunctions[7] = -32.0 * r * s + 160.0 * r * r * s + 32.0 * r * s * s - 128.0 * r * r * s * s - 128.0 * r * r * r * s;
    shapeFunctions[8] = +5.33333333333 * r * s - 32.0 * r * r * s + 42.6666666667 * r * r * r * s;
    shapeFunctions[9] = -12.0 * s + 28.0 * r * s + 76.0 * s * s - 16.0 * r * r * s - 144.0 * r * s * s - 128.0 * s * s * s + 64.0 * s * s * s * s + 128.0 * r * s * s * s + 64.0 * r * r * s * s;
    shapeFunctions[10] = -32.0 * r * s + 32.0 * r * r * s + 160.0 * r * s * s - 128.0 * r * s * s * s - 128.0 * r * r * s * s;
    shapeFunctions[11] = +4.0 * r * s - 16.0 * r * r * s - 16.0 * r * s * s + 64.0 * r * r * s * s;
    shapeFunctions[12] = +5.33333333333 * s - 5.33333333333 * r * s - 37.3333333333 * s * s + 32.0 * r * s * s + 74.6666666667 * s * s * s - 42.6666666667 * s * s * s * s - 42.6666666667 * r * s * s * s;
    shapeFunctions[13] = +5.33333333333 * r * s - 32.0 * r * s * s + 42.6666666667 * r * s * s * s;
    shapeFunctions[14] = -1.0 * s + 7.33333333333 * s * s - 16.0 * s * s * s + 10.6666666667 * s * s * s * s;

    return shapeFunctions;
}

Eigen::Matrix<double, 15, 2> DerivativeShapeFunctionsTriangleOrder4(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 15, 2> derivativeShapeFunctions;
    double r(rCoordinates(0));
    double s(rCoordinates(1));

    derivativeShapeFunctions(0, 0) = -8.33333333333 + 46.6666666667 * r + 46.6666666667 * s - 80.0 * r * r - 160.0 * r * s - 80.0 * s * s + 42.6666666667 * s * s * s + 128.0 * r * s * s + 128.0 * r * r * s + 42.6666666667 * r * r * r;
    derivativeShapeFunctions(0, 1) = -8.33333333333 + 46.6666666667 * r + 46.6666666667 * s - 80.0 * r * r - 160.0 * r * s - 80.0 * s * s + 42.6666666667 * s * s * s + 128.0 * r * s * s + 128.0 * r * r * s + 42.6666666667 * r * r * r;

    derivativeShapeFunctions(1, 0) = +16.0 - 138.666666667 * r - 69.3333333333 * s + 288.0 * r * r + 384.0 * r * s + 96.0 * s * s - 42.6666666667 * s * s * s - 256.0 * r * s * s - 384.0 * r * r * s - 170.666666667 * r * r * r;
    derivativeShapeFunctions(1, 1) = -69.3333333333 * r + 192.0 * r * r + 192.0 * r * s - 128.0 * r * s * s - 256.0 * r * r * s - 128.0 * r * r * r;

    derivativeShapeFunctions(2, 0) = -12.0 + 152.0 * r + 28.0 * s - 384.0 * r * r - 288.0 * r * s - 16.0 * s * s + 128.0 * r * s * s + 384.0 * r * r * s + 256.0 * r * r * r;
    derivativeShapeFunctions(2, 1) = +28.0 * r - 144.0 * r * r - 32.0 * r * s + 128.0 * r * r * s + 128.0 * r * r * r;

    derivativeShapeFunctions(3, 0) = +5.33333333333 - 74.6666666667 * r - 5.33333333333 * s + 224.0 * r * r + 64.0 * r * s - 128.0 * r * r * s - 170.666666667 * r * r * r;
    derivativeShapeFunctions(3, 1) = -5.33333333333 * r + 32.0 * r * r - 42.6666666667 * r * r * r;

    derivativeShapeFunctions(4, 0) = -1.0 + 14.6666666667 * r - 48.0 * r * r + 42.6666666667 * r * r * r;
    derivativeShapeFunctions(4, 1) = 0.;

    derivativeShapeFunctions(5, 0) = -69.3333333333 * s + 192.0 * r * s + 192.0 * s * s - 128.0 * s * s * s - 256.0 * r * s * s - 128.0 * r * r * s;
    derivativeShapeFunctions(5, 1) = +16.0 - 69.3333333333 * r - 138.666666667 * s + 96.0 * r * r + 384.0 * r * s + 288.0 * s * s - 170.666666667 * s * s * s - 384.0 * r * s * s - 256.0 * r * r * s - 42.6666666667 * r * r * r;

    derivativeShapeFunctions(6, 0) = +96.0 * s - 448.0 * r * s - 224.0 * s * s + 128.0 * s * s * s + 512.0 * r * s * s + 384.0 * r * r * s;
    derivativeShapeFunctions(6, 1) = +96.0 * r - 224.0 * r * r - 448.0 * r * s + 384.0 * r * s * s + 512.0 * r * r * s + 128.0 * r * r * r;

    derivativeShapeFunctions(7, 0) = -32.0 * s + 320.0 * r * s + 32.0 * s * s - 256.0 * r * s * s - 384.0 * r * r * s;
    derivativeShapeFunctions(7, 1) = -32.0 * r + 160.0 * r * r + 64.0 * r * s - 256.0 * r * r * s - 128.0 * r * r * r;

    derivativeShapeFunctions(8, 0) = +5.33333333333 * s - 64.0 * r * s + 128.0 * r * r * s;
    derivativeShapeFunctions(8, 1) = +5.33333333333 * r - 32.0 * r * r + 42.6666666667 * r * r * r;

    derivativeShapeFunctions(9, 0) = +28.0 * s - 32.0 * r * s - 144.0 * s * s + 128.0 * s * s * s + 128.0 * r * s * s;
    derivativeShapeFunctions(9, 1) = -12.0 + 28.0 * r + 152.0 * s - 16.0 * r * r - 288.0 * r * s - 384.0 * s * s + 256.0 * s * s * s + 384.0 * r * s * s + 128.0 * r * r * s;

    derivativeShapeFunctions(10, 0) = -32.0 * s + 64.0 * r * s + 160.0 * s * s - 128.0 * s * s * s - 256.0 * r * s * s;
    derivativeShapeFunctions(10, 1) = -32.0 * r + 32.0 * r * r + 320.0 * r * s - 384.0 * r * s * s - 256.0 * r * r * s;

    derivativeShapeFunctions(11, 0) = +4.0 * s - 32.0 * r * s - 16.0 * s * s + 128.0 * r * s * s;
    derivativeShapeFunctions(11, 1) = +4.0 * r - 16.0 * r * r - 32.0 * r * s + 128.0 * r * r * s;

    derivativeShapeFunctions(12, 0) = -5.33333333333 * s + 32.0 * s * s - 42.6666666667 * s * s * s;
    derivativeShapeFunctions(12, 1) = +5.33333333333 - 5.33333333333 * r - 74.6666666667 * s + 64.0 * r * s + 224.0 * s * s - 170.666666667 * s * s * s - 128.0 * r * s * s;

    derivativeShapeFunctions(13, 0) = +5.33333333333 * s - 32.0 * s * s + 42.6666666667 * s * s * s;
    derivativeShapeFunctions(13, 1) = +5.33333333333 * r - 64.0 * r * s + 128.0 * r * s * s;

    derivativeShapeFunctions(14, 0) = 0.;
    derivativeShapeFunctions(14, 1) = -1.0 + 14.6666666667 * s - 48.0 * s * s + 42.6666666667 * s * s * s;

    return derivativeShapeFunctions;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 2, 1> NodeCoordinatesQuadOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Vector2d(-1.0, -1.0);
    case 1: return Eigen::Vector2d(1.0, -1.0);
    case 2: return Eigen::Vector2d(1.0, 1.0);
    case 3: return Eigen::Vector2d(-1.0, 1.0);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..3)");
    }
}

Eigen::Matrix<double, 4, 1> ShapeFunctionsQuadOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> shapeFunctions;
    shapeFunctions[0] = 0.25 * (1. - rCoordinates(0)) * (1. - rCoordinates(1));
    shapeFunctions[1] = 0.25 * (1. + rCoordinates(0)) * (1. - rCoordinates(1));
    shapeFunctions[2] = 0.25 * (1. + rCoordinates(0)) * (1. + rCoordinates(1));
    shapeFunctions[3] = 0.25 * (1. - rCoordinates(0)) * (1. + rCoordinates(1));
    return shapeFunctions;
}

Eigen::Matrix<double, 4, 2> DerivativeShapeFunctionsQuadOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 2> derivativeShapeFunctions;
    derivativeShapeFunctions(0, 0) = -0.25 * (1. - rCoordinates(1));
    derivativeShapeFunctions(0, 1) = -0.25 * (1. - rCoordinates(0));

    derivativeShapeFunctions(1, 0) = +0.25 * (1. - rCoordinates(1));
    derivativeShapeFunctions(1, 1) = -0.25 * (1. + rCoordinates(0));

    derivativeShapeFunctions(2, 0) = +0.25 * (1. + rCoordinates(1));
    derivativeShapeFunctions(2, 1) = +0.25 * (1. + rCoordinates(0));

    derivativeShapeFunctions(3, 0) = -0.25 * (1. + rCoordinates(1));
    derivativeShapeFunctions(3, 1) = +0.25 * (1. - rCoordinates(0));

    return derivativeShapeFunctions;
}

Eigen::Matrix<double, 2, 1> NodeCoordinatesQuadOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Vector2d(-1.0, -1.0);
    case 1: return Eigen::Vector2d(1.0, -1.0);
    case 2: return Eigen::Vector2d(1.0, 1.0);
    case 3: return Eigen::Vector2d(-1.0, 1.0);
    case 4: return Eigen::Vector2d(0.0, -1.0);
    case 5: return Eigen::Vector2d(1.0, 0.0);
    case 6: return Eigen::Vector2d(0.0, 1.0);
    case 7: return Eigen::Vector2d(-1.0, 0.0);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..7)");
    }
}

Eigen::Matrix<double, 8, 1> ShapeFunctionsQuadOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 8, 1> shapeFunctions;
    double r = rCoordinates(0);
    double s = rCoordinates(1);

    shapeFunctions[0] = -0.25 * (r - 1) * (s - 1) * (r + s + 1);
    shapeFunctions[1] = 0.25 * (r + 1) * (s - 1) * (-r + s + 1);
    shapeFunctions[2] = 0.25 * (r + 1) * (s + 1) * (r + s - 1);
    shapeFunctions[3] = 0.25 * (r - 1) * (s + 1) * (r - s + 1);
    shapeFunctions[4] = 0.5 * (r * r - 1) * (s - 1);
    shapeFunctions[5] = -0.5 * (r + 1) * (s * s - 1);
    shapeFunctions[6] = -0.5 * (r * r - 1) * (s + 1);
    shapeFunctions[7] = 0.5 * (r - 1) * (s * s - 1);

    return shapeFunctions;
}

Eigen::Matrix<double, 8, 2> DerivativeShapeFunctionsQuadOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 8, 2> derivativeShapeFunctions;
    double r = rCoordinates(0);
    double s = rCoordinates(1);

    derivativeShapeFunctions(0, 0) = (2 * r + s) * (-0.25 * s + 0.25);
    derivativeShapeFunctions(0, 1) = (-0.25 * r + 0.25) * (r + 2 * s);

    derivativeShapeFunctions(1, 0) = (-0.5 * r + 0.25 * s) * (s - 1);
    derivativeShapeFunctions(1, 1) = 0.25 * (-r + 2 * s) * (r + 1);

    derivativeShapeFunctions(2, 0) = 0.25 * (2 * r + s) * (s + 1);
    derivativeShapeFunctions(2, 1) = 0.25 * (r + 1) * (r + 2 * s);

    derivativeShapeFunctions(3, 0) = 0.25 * (2 * r - s) * (s + 1);
    derivativeShapeFunctions(3, 1) = (0.25 * r - 0.5 * s) * (r - 1);

    derivativeShapeFunctions(4, 0) = r * (s - 1);
    derivativeShapeFunctions(4, 1) = 0.5 * r * r - 0.5;

    derivativeShapeFunctions(5, 0) = -0.5 * s * s + 0.5;
    derivativeShapeFunctions(5, 1) = -s * (r + 1);

    derivativeShapeFunctions(6, 0) = -r * (s + 1);
    derivativeShapeFunctions(6, 1) = -0.5 * r * r + 0.5;

    derivativeShapeFunctions(7, 0) = 0.5 * s * s - 0.5;
    derivativeShapeFunctions(7, 1) = s * (r - 1);

    return derivativeShapeFunctions;
}

Eigen::Matrix<double, 2, 1> NodeCoordinatesQuadSpectralOrder2(int rNodeIndex)
{
    const int d = 3;

    assert(rNodeIndex >= 0);
    assert(rNodeIndex < 9);

    double cX = ShapeFunctions1D::NodeCoordinatesTrussOrder2(rNodeIndex % d)(0, 0);
    double cY = ShapeFunctions1D::NodeCoordinatesTrussOrder2(rNodeIndex / d)(0, 0);

    return Eigen::Vector2d({cX, cY});
}

Eigen::Matrix<double, 9, 1> ShapeFunctionsQuadSpectralOrder2(const Eigen::VectorXd& rCoordinates)
{
    const int D = 3;
    const int DxD = 9;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussOrder2(cX);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussOrder2(cY);

    Eigen::Matrix<double, D, D> shapeFunctionsMatrix = shapeFunctions1Dx * shapeFunctions1Dy.transpose();
    return Eigen::Map<Eigen::Matrix<double, DxD, 1>>(shapeFunctionsMatrix.data(), DxD);
}

Eigen::Matrix<double, 9, 2> DerivativeShapeFunctionsQuadSpectralOrder2(const Eigen::VectorXd& rCoordinates)
{
    const int D = 3;
    const int DxD = 9;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussOrder2(cX);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussOrder2(cY);
    const Eigen::Matrix<double, D, 1>& derShapeFunctions1Dx = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(cX);
    const Eigen::Matrix<double, D, 1>& derShapeFunctions1Dy = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(cY);

    Eigen::Matrix<double, 9, 2> derivativeShapeFunctions;

    Eigen::Matrix<double, D, D> dNdXi = derShapeFunctions1Dx * shapeFunctions1Dy.transpose();
    Eigen::Matrix<double, D, D> dNdEta = shapeFunctions1Dx * derShapeFunctions1Dy.transpose();

    derivativeShapeFunctions.block<DxD, 1>(0, 0) = Eigen::Map<Eigen::Matrix<double, DxD, 1>>(dNdXi.data(), DxD);
    derivativeShapeFunctions.block<DxD, 1>(0, 1) = Eigen::Map<Eigen::Matrix<double, DxD, 1>>(dNdEta.data(), DxD);

    //    int theDerShapeFunction(0);
//    for (int county = 0; county < D; county ++)
//    {
//        for (int countx = 0; countx < D; countx ++)
//        {
//            derivativeShapeFunctions(theDerShapeFunction,0) = derShapeFunctions1Dx[countx]*shapeFunctions1Dy[county];
//            derivativeShapeFunctions(theDerShapeFunction,1) = shapeFunctions1Dx[countx]*derShapeFunctions1Dy[county];
//            theDerShapeFunction++;
//        }
//    }

    return derivativeShapeFunctions;

}
Eigen::Matrix<double, 2, 1> NodeCoordinatesQuadSpectralOrder3(int rNodeIndex)
{
    const int d = 4;

    assert(rNodeIndex >= 0);
    assert(rNodeIndex < 16);

    double cX = ShapeFunctions1D::NodeCoordinatesTrussSpectralOrder3(rNodeIndex % d)(0, 0);
    double cY = ShapeFunctions1D::NodeCoordinatesTrussSpectralOrder3(rNodeIndex / d)(0, 0);

    return Eigen::Vector2d({cX, cY});
}

Eigen::Matrix<double, 16, 1> ShapeFunctionsQuadSpectralOrder3(const Eigen::VectorXd& rCoordinates)
{
    const int D = 4;
    const int DxD = 16;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(cX);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(cY);

    Eigen::Matrix<double, D, D> shapeFunctionsMatrix = shapeFunctions1Dx * shapeFunctions1Dy.transpose();
    return Eigen::Map<Eigen::Matrix<double, DxD, 1>>(shapeFunctionsMatrix.data(), DxD);

}

Eigen::Matrix<double, 16, 2> DerivativeShapeFunctionsQuadSpectralOrder3(const Eigen::VectorXd& rCoordinates)
{
    const int D = 4;
    const int DxD = 16;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(cX);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(cY);
    const Eigen::Matrix<double, D, 1>& derShapeFunctions1Dx = ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder3(cX);
    const Eigen::Matrix<double, D, 1>& derShapeFunctions1Dy = ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder3(cY);

    Eigen::Matrix<double, DxD, 2> derivativeShapeFunctions;

    Eigen::Matrix<double, D, D> dNdXi = derShapeFunctions1Dx * shapeFunctions1Dy.transpose();
    Eigen::Matrix<double, D, D> dNdEta = shapeFunctions1Dx * derShapeFunctions1Dy.transpose();

    derivativeShapeFunctions.block<DxD, 1>(0, 0) = Eigen::Map<Eigen::Matrix<double, DxD, 1>>(dNdXi.data(), DxD);
    derivativeShapeFunctions.block<DxD, 1>(0, 1) = Eigen::Map<Eigen::Matrix<double, DxD, 1>>(dNdEta.data(), DxD);

    //    int theDerShapeFunction(0);
//    for (int county = 0; county < D; county ++)
//    {
//        for (int countx = 0; countx < D; countx ++)
//        {
//            derivativeShapeFunctions(theDerShapeFunction,0) = derShapeFunctions1Dx[countx]*shapeFunctions1Dy[county];
//            derivativeShapeFunctions(theDerShapeFunction,1) = shapeFunctions1Dx[countx]*derShapeFunctions1Dy[county];
//            theDerShapeFunction++;
//        }
//    }

    return derivativeShapeFunctions;
}

Eigen::Matrix<double, 2, 1> NodeCoordinatesQuadSpectralOrder4(int rNodeIndex)
{
    const int d = 5;

    assert(rNodeIndex >= 0);
    assert(rNodeIndex < 25);

    double cX = ShapeFunctions1D::NodeCoordinatesTrussSpectralOrder4(rNodeIndex % d)(0, 0);
    double cY = ShapeFunctions1D::NodeCoordinatesTrussSpectralOrder4(rNodeIndex / d)(0, 0);

    return Eigen::Vector2d({cX, cY});
}

Eigen::Matrix<double, 25, 1> ShapeFunctionsQuadSpectralOrder4(const Eigen::VectorXd& rCoordinates)
{
    const int D = 5;
    const int DxD = 25;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(cX);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(cY);

    Eigen::Matrix<double, D, D> shapeFunctionsMatrix = shapeFunctions1Dx * shapeFunctions1Dy.transpose();
    return Eigen::Map<Eigen::Matrix<double, DxD, 1>>(shapeFunctionsMatrix.data(), DxD);

}

Eigen::Matrix<double, 25, 2> DerivativeShapeFunctionsQuadSpectralOrder4(const Eigen::VectorXd& rCoordinates)
{
    const int D = 5;
    const int DxD = 25;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(cX);
    const Eigen::Matrix<double, D, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(cY);
    const Eigen::Matrix<double, D, 1>& derShapeFunctions1Dx = ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder4(cX);
    const Eigen::Matrix<double, D, 1>& derShapeFunctions1Dy = ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder4(cY);

    Eigen::Matrix<double, DxD, 2> derivativeShapeFunctions;

    Eigen::Matrix<double, D, D> dNdXi = derShapeFunctions1Dx * shapeFunctions1Dy.transpose();
    Eigen::Matrix<double, D, D> dNdEta = shapeFunctions1Dx * derShapeFunctions1Dy.transpose();

    derivativeShapeFunctions.block<DxD, 1>(0, 0) = Eigen::Map<Eigen::Matrix<double, DxD, 1>>(dNdXi.data(), DxD);
    derivativeShapeFunctions.block<DxD, 1>(0, 1) = Eigen::Map<Eigen::Matrix<double, DxD, 1>>(dNdEta.data(), DxD);

    //    int theDerShapeFunction(0);
//    for (int county = 0; county < D; county ++)
//    {
//        for (int countx = 0; countx < D; countx ++)
//        {
//            derivativeShapeFunctions(theDerShapeFunction,0) = derShapeFunctions1Dx[countx]*shapeFunctions1Dy[county];
//            derivativeShapeFunctions(theDerShapeFunction,1) = shapeFunctions1Dx[countx]*derShapeFunctions1Dy[county];
//            theDerShapeFunction++;
//        }
//    }

    return derivativeShapeFunctions;
}

}

namespace ShapeFunctions3D
{

Eigen::Matrix<double, 3, 1> NodeCoordinatesTetrahedronOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Vector3d(0.0, 0.0, 0.0);
    case 1: return Eigen::Vector3d(1.0, 0.0, 0.0);
    case 2: return Eigen::Vector3d(0.0, 1.0, 0.0);
    case 3: return Eigen::Vector3d(0.0, 0.0, 1.0);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..4)");
    }
}

Eigen::Matrix<double, 4, 1> ShapeFunctionsTetrahedronOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 1> shapeFunctions;
    shapeFunctions[0] = 1. - rCoordinates.sum();
    shapeFunctions[1] = rCoordinates(0);
    shapeFunctions[2] = rCoordinates(1);
    shapeFunctions[3] = rCoordinates(2);
    return shapeFunctions;
}

Eigen::Matrix<double, 4, 3> DerivativeShapeFunctionsTetrahedronOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 4, 3> derivativeShapeFunctions;
    derivativeShapeFunctions(0, 0) = -1;
    derivativeShapeFunctions(0, 1) = -1;
    derivativeShapeFunctions(0, 2) = -1;

    derivativeShapeFunctions(1, 0) = 1;
    derivativeShapeFunctions(1, 1) = 0;
    derivativeShapeFunctions(1, 2) = 0;

    derivativeShapeFunctions(2, 0) = 0;
    derivativeShapeFunctions(2, 1) = 1;
    derivativeShapeFunctions(2, 2) = 0;

    derivativeShapeFunctions(3, 0) = 0;
    derivativeShapeFunctions(3, 1) = 0;
    derivativeShapeFunctions(3, 2) = 1;

    return derivativeShapeFunctions;
}

Eigen::Matrix<double, 3, 1> NodeCoordinatesTetrahedronOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {

    case 0: return Eigen::Vector3d(0.0, 0.0, 0.0);
    case 1: return Eigen::Vector3d(1.0, 0.0, 0.0);
    case 2: return Eigen::Vector3d(0.0, 1.0, 0.0);
    case 3: return Eigen::Vector3d(0.0, 0.0, 1.0);
    case 4: return Eigen::Vector3d(0.5, 0.0, 0.0);
    case 5: return Eigen::Vector3d(0.5, 0.5, 0.0);
    case 6: return Eigen::Vector3d(0.0, 0.5, 0.0);
    case 7: return Eigen::Vector3d(0.0, 0.0, 0.5);
    case 8: return Eigen::Vector3d(0.0, 0.5, 0.5);
    case 9: return Eigen::Vector3d(0.5, 0.0, 0.5);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..9)");
    }

}

Eigen::Matrix<double, 10, 1> ShapeFunctionsTetrahedronOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 10, 1> shapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];
    double t = rCoordinates[2];

    // node 0 (0,0,0)
    shapeFunctions[0] = 1. - 3. * (r + s + t) + 2. * (r * r + s * s + t * t) + 4. * (r * s + r * t + s * t);

    // node 1 (1,0,0)
    shapeFunctions[1] = -r + 2. * r * r;

    // node 2 (0,1,0)
    shapeFunctions[2] = -s + 2. * s * s;

    // node 3 (0,0,1)
    shapeFunctions[3] = -t + 2. * t * t;

    // node 4 (0.5,0,0)
    shapeFunctions[4] = 4. * r * (1. - r - s - t);

    // node 5 (0.5,0.5,0)
    shapeFunctions[5] = 4. * r * s;

    // node 6 (0,0.5,0)
    shapeFunctions[6] = 4. * s * (1. - r - s - t);

    // node 7 (0,0,0.5)
    shapeFunctions[7] = 4. * t * (1. - r - s - t);

    // node 8 (0,0.5,0.5)
    shapeFunctions[8] = 4. * s * t;

    // node 9 (0.5,0,0.5)
    shapeFunctions[9] = 4. * r * t;

    return shapeFunctions;
}

Eigen::Matrix<double, 10, 3> DerivativeShapeFunctionsTetrahedronOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 10, 3> derivativeShapeFunctions;
    double r = rCoordinates[0];
    double s = rCoordinates[1];
    double t = rCoordinates[2];

    derivativeShapeFunctions(0, 0) = -3. + 4. * (r + s + t);
    derivativeShapeFunctions(0, 1) = -3. + 4. * (r + s + t);
    derivativeShapeFunctions(0, 2) = -3. + 4. * (r + s + t);

    derivativeShapeFunctions(1, 0) = -1. + 4. * r;
    derivativeShapeFunctions(1, 1) = 0;
    derivativeShapeFunctions(1, 2) = 0;

    derivativeShapeFunctions(2, 0) = 0;
    derivativeShapeFunctions(2, 1) = -1. + 4. * s;
    derivativeShapeFunctions(2, 2) = 0;

    derivativeShapeFunctions(3, 0) = 0;
    derivativeShapeFunctions(3, 1) = 0;
    derivativeShapeFunctions(3, 2) = -1. + 4. * t;

    derivativeShapeFunctions(4, 0) = 4. - 8. * r - 4. * s - 4. * t;
    derivativeShapeFunctions(4, 1) = -4. * r;
    derivativeShapeFunctions(4, 2) = -4. * r;

    derivativeShapeFunctions(5, 0) = 4. * s;
    derivativeShapeFunctions(5, 1) = 4. * r;
    derivativeShapeFunctions(5, 2) = 0.;

    derivativeShapeFunctions(6, 0) = -4. * s;
    derivativeShapeFunctions(6, 1) = 4. - 8. * s - 4. * r - 4. * t;
    derivativeShapeFunctions(6, 2) = -4. * s;

    derivativeShapeFunctions(7, 0) = -4. * t;
    derivativeShapeFunctions(7, 1) = -4. * t;
    derivativeShapeFunctions(7, 2) = 4. - 8. * t - 4. * r - 4. * s;

    derivativeShapeFunctions(8, 0) = 0.;
    derivativeShapeFunctions(8, 1) = 4. * t;
    derivativeShapeFunctions(8, 2) = 4. * s;

    derivativeShapeFunctions(9, 0) = 4. * t;
    derivativeShapeFunctions(9, 1) = 0;
    derivativeShapeFunctions(9, 2) = 4. * r;

    return derivativeShapeFunctions;
}

Eigen::Matrix<double, 3, 1> NodeCoordinatesBrickOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {

    case 0: return Eigen::Vector3d(-1., -1., -1.);
    case 1: return Eigen::Vector3d(1., -1., -1.);
    case 2: return Eigen::Vector3d(1., 1., -1.);
    case 3: return Eigen::Vector3d(-1., 1., -1.);
    case 4: return Eigen::Vector3d(-1., -1., 1.);
    case 5: return Eigen::Vector3d(1., -1., 1.);
    case 6: return Eigen::Vector3d(1., 1., 1.);
    case 7: return Eigen::Vector3d(-1., 1., 1.);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..7)");
    }
}

Eigen::Matrix<double, 8, 1> ShapeFunctionsBrickOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 8, 1> shapeFunctions;

    double plus_r = 1.0 + rCoordinates[0];
    double plus_s = 1.0 + rCoordinates[1];
    double plus_t = 1.0 + rCoordinates[2];

    double minus_r = 1.0 - rCoordinates[0];
    double minus_s = 1.0 - rCoordinates[1];
    double minus_t = 1.0 - rCoordinates[2];

    shapeFunctions[0] = 0.125 * minus_r * minus_s * minus_t;
    shapeFunctions[1] = 0.125 * plus_r * minus_s * minus_t;
    shapeFunctions[2] = 0.125 * plus_r * plus_s * minus_t;
    shapeFunctions[3] = 0.125 * minus_r * plus_s * minus_t;
    shapeFunctions[4] = 0.125 * minus_r * minus_s * plus_t;
    shapeFunctions[5] = 0.125 * plus_r * minus_s * plus_t;
    shapeFunctions[6] = 0.125 * plus_r * plus_s * plus_t;
    shapeFunctions[7] = 0.125 * minus_r * plus_s * plus_t;

    return shapeFunctions;
}

Eigen::Matrix<double, 8, 3> DerivativeShapeFunctionsBrickOrder1(const Eigen::VectorXd& rCoordinates)
{
    double plus_r = 1.0 + rCoordinates[0];
    double plus_s = 1.0 + rCoordinates[1];
    double plus_t = 1.0 + rCoordinates[2];

    double minus_r = 1.0 - rCoordinates[0];
    double minus_s = 1.0 - rCoordinates[1];
    double minus_t = 1.0 - rCoordinates[2];

    Eigen::Matrix<double, 8, 3> derivativeShapeFunctions;

    derivativeShapeFunctions(0, 0) = -0.125 * minus_s * minus_t;
    derivativeShapeFunctions(0, 1) = -0.125 * minus_r * minus_t;
    derivativeShapeFunctions(0, 2) = -0.125 * minus_r * minus_s;

    derivativeShapeFunctions(1, 0) = 0.125 * minus_s * minus_t;
    derivativeShapeFunctions(1, 1) = -0.125 * plus_r * minus_t;
    derivativeShapeFunctions(1, 2) = -0.125 * plus_r * minus_s;

    derivativeShapeFunctions(2, 0) = 0.125 * plus_s * minus_t;
    derivativeShapeFunctions(2, 1) = 0.125 * plus_r * minus_t;
    derivativeShapeFunctions(2, 2) = -0.125 * plus_r * plus_s;

    derivativeShapeFunctions(3, 0) = -0.125 * plus_s * minus_t;
    derivativeShapeFunctions(3, 1) = 0.125 * minus_r * minus_t;
    derivativeShapeFunctions(3, 2) = -0.125 * minus_r * plus_s;

    derivativeShapeFunctions(4, 0) = -0.125 * minus_s * plus_t;
    derivativeShapeFunctions(4, 1) = -0.125 * minus_r * plus_t;
    derivativeShapeFunctions(4, 2) = 0.125 * minus_r * minus_s;

    derivativeShapeFunctions(5, 0) = 0.125 * minus_s * plus_t;
    derivativeShapeFunctions(5, 1) = -0.125 * plus_r * plus_t;
    derivativeShapeFunctions(5, 2) = 0.125 * plus_r * minus_s;

    derivativeShapeFunctions(6, 0) = 0.125 * plus_s * plus_t;
    derivativeShapeFunctions(6, 1) = 0.125 * plus_r * plus_t;
    derivativeShapeFunctions(6, 2) = 0.125 * plus_r * plus_s;

    derivativeShapeFunctions(7, 0) = -0.125 * plus_s * plus_t;
    derivativeShapeFunctions(7, 1) = 0.125 * minus_r * plus_t;
    derivativeShapeFunctions(7, 2) = 0.125 * minus_r * plus_s;

    return derivativeShapeFunctions;
}

Eigen::Matrix<double, 3, 1> NodeCoordinatesBrickOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {

    case 0: return Eigen::Vector3d(-1., -1., -1.);
    case 1: return Eigen::Vector3d(1., -1., -1.);
    case 2: return Eigen::Vector3d(1., 1., -1.);
    case 3: return Eigen::Vector3d(-1., 1., -1.);
    case 4: return Eigen::Vector3d(-1., -1., 1.);
    case 5: return Eigen::Vector3d(1., -1., 1.);
    case 6: return Eigen::Vector3d(1., 1., 1.);
    case 7: return Eigen::Vector3d(-1., 1., 1.);

    case 8: return Eigen::Vector3d(0., -1., -1.);
    case 9: return Eigen::Vector3d(1., 0., -1.);
    case 10: return Eigen::Vector3d(0., 1., -1.);
    case 11: return Eigen::Vector3d(-1., 0., -1.);

    case 12: return Eigen::Vector3d(-1., -1., 0.);
    case 13: return Eigen::Vector3d(1., -1., 0.);
    case 14: return Eigen::Vector3d(1., 1., 0.);
    case 15: return Eigen::Vector3d(-1., 1., 0.);

    case 16: return Eigen::Vector3d(0., -1., 1.);
    case 17: return Eigen::Vector3d(1., 0., 1.);
    case 18: return Eigen::Vector3d(0., 1., 1.);
    case 19: return Eigen::Vector3d(-1., 0., 1.);

    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..19)");
    }
}

Eigen::Matrix<double, 20, 1> ShapeFunctionsBrickOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 20, 1> shapeFunctions;

    double r = rCoordinates[0];
    double s = rCoordinates[1];
    double t = rCoordinates[2];

    double plus_r = 1.0 + r;
    double plus_s = 1.0 + s;
    double plus_t = 1.0 + t;

    double minus_r = 1.0 - r;
    double minus_s = 1.0 - s;
    double minus_t = 1.0 - t;

    double mid_r = 1.0 - r * r;
    double mid_s = 1.0 - s * s;
    double mid_t = 1.0 - t * t;

    shapeFunctions[0] = 0.125 * minus_r * minus_s * minus_t * (-r - s - t - 2);
    shapeFunctions[1] = 0.125 * plus_r * minus_s * minus_t * (r - s - t - 2);
    shapeFunctions[2] = 0.125 * plus_r * plus_s * minus_t * (r + s - t - 2);
    shapeFunctions[3] = 0.125 * minus_r * plus_s * minus_t * (-r + s - t - 2);
    shapeFunctions[4] = 0.125 * minus_r * minus_s * plus_t * (-r - s + t - 2);
    shapeFunctions[5] = 0.125 * plus_r * minus_s * plus_t * (r - s + t - 2);
    shapeFunctions[6] = 0.125 * plus_r * plus_s * plus_t * (r + s + t - 2);
    shapeFunctions[7] = 0.125 * minus_r * plus_s * plus_t * (-r + s + t - 2);

    shapeFunctions[8] = 0.25 * mid_r * minus_s * minus_t;
    shapeFunctions[9] = 0.25 * plus_r * mid_s * minus_t;
    shapeFunctions[10] = 0.25 * mid_r * plus_s * minus_t;
    shapeFunctions[11] = 0.25 * minus_r * mid_s * minus_t;

    shapeFunctions[12] = 0.25 * minus_r * minus_s * mid_t;
    shapeFunctions[13] = 0.25 * plus_r * minus_s * mid_t;
    shapeFunctions[14] = 0.25 * plus_r * plus_s * mid_t;
    shapeFunctions[15] = 0.25 * minus_r * plus_s * mid_t;

    shapeFunctions[16] = 0.25 * mid_r * minus_s * plus_t;
    shapeFunctions[17] = 0.25 * plus_r * mid_s * plus_t;
    shapeFunctions[18] = 0.25 * mid_r * plus_s * plus_t;
    shapeFunctions[19] = 0.25 * minus_r * mid_s * plus_t;

    return shapeFunctions;
}

Eigen::Matrix<double, 20, 3> DerivativeShapeFunctionsBrickOrder2(const Eigen::VectorXd& rCoordinates)
{

    double r = rCoordinates[0];
    double s = rCoordinates[1];
    double t = rCoordinates[2];

    double plus_r = 1.0 + r;
    double plus_s = 1.0 + s;
    double plus_t = 1.0 + t;

    double minus_r = 1.0 - r;
    double minus_s = 1.0 - s;
    double minus_t = 1.0 - t;

    double mid_r = 1.0 - r * r;
    double mid_s = 1.0 - s * s;
    double mid_t = 1.0 - t * t;

    Eigen::Matrix<double, 20, 3> derivativeShapeFunctions;

//    shapeFunctions[0] = 0.125 * minus_r * minus_s * minus_t * (-r-s-t-2);
    derivativeShapeFunctions(0, 0) = 0.125 * minus_s * minus_t * (2 * r + s + t + 1);
    derivativeShapeFunctions(0, 1) = 0.125 * minus_r * minus_t * (r + 2 * s + t + 1);
    derivativeShapeFunctions(0, 2) = 0.125 * minus_r * minus_s * (r + s + 2 * t + 1);

//    shapeFunctions[1] = 0.125 * plus_r  * minus_s * minus_t * ( r-s-t-2);
    derivativeShapeFunctions(1, 0) = 0.125 * minus_s * minus_t * (2 * r - s - t - 1);
    derivativeShapeFunctions(1, 1) = -0.125 * plus_r * minus_t * (r - 2 * s - t - 1);
    derivativeShapeFunctions(1, 2) = -0.125 * plus_r * minus_s * (r - s - 2 * t - 1);

//    shapeFunctions[2] = 0.125 * plus_r  * plus_s  * minus_t * ( r+s-t-2);
    derivativeShapeFunctions(2, 0) = 0.125 * plus_s * minus_t * (2 * r + s - t - 1);
    derivativeShapeFunctions(2, 1) = 0.125 * plus_r * minus_t * (r + 2 * s - t - 1);
    derivativeShapeFunctions(2, 2) = -0.125 * plus_r * plus_s * (r + s - 2 * t - 1);

//    shapeFunctions[3] = 0.125 * minus_r * plus_s  * minus_t * (-r+s-t-2);
    derivativeShapeFunctions(3, 0) = 0.125 * plus_s * minus_t * (2 * r - s + t + 1);
    derivativeShapeFunctions(3, 1) = -0.125 * minus_r * minus_t * (r - 2 * s + t + 1);
    derivativeShapeFunctions(3, 2) = 0.125 * minus_r * plus_s * (r - s + 2 * t + 1);

//    shapeFunctions[4] = 0.125 * minus_r * minus_s *  plus_t * (-r-s+t-2);
    derivativeShapeFunctions(4, 0) = 0.125 * minus_s * plus_t * (2 * r + s - t + 1);
    derivativeShapeFunctions(4, 1) = 0.125 * minus_r * plus_t * (r + 2 * s - t + 1);
    derivativeShapeFunctions(4, 2) = -0.125 * minus_r * minus_s * (r + s - 2 * t + 1);

//    shapeFunctions[5] = 0.125 * plus_r  * minus_s *  plus_t * ( r-s+t-2);
    derivativeShapeFunctions(5, 0) = 0.125 * minus_s * plus_t * (2 * r - s + t - 1);
    derivativeShapeFunctions(5, 1) = -0.125 * plus_r * plus_t * (r - 2 * s + t - 1);
    derivativeShapeFunctions(5, 2) = 0.125 * plus_r * minus_s * (r - s + 2 * t - 1);

//    shapeFunctions[6] = 0.125 * plus_r  * plus_s  *  plus_t * ( r+s+t-2);
    derivativeShapeFunctions(6, 0) = 0.125 * plus_s * plus_t * (2 * r + s + t - 1);
    derivativeShapeFunctions(6, 1) = 0.125 * plus_r * plus_t * (r + 2 * s + t - 1);
    derivativeShapeFunctions(6, 2) = 0.125 * plus_r * plus_s * (r + s + 2 * t - 1);

//    shapeFunctions[7] = 0.125 * minus_r * plus_s  *  plus_t * (-r+s+t-2);
    derivativeShapeFunctions(7, 0) = 0.125 * plus_s * plus_t * (2 * r - s - t + 1);
    derivativeShapeFunctions(7, 1) = -0.125 * minus_r * plus_t * (r - 2 * s - t + 1);
    derivativeShapeFunctions(7, 2) = -0.125 * minus_r * plus_s * (r - s - 2 * t + 1);

//    shapeFunctions[8] =  0.25 *   mid_r * minus_s * minus_t;
    derivativeShapeFunctions(8, 0) = -0.5 * r * minus_s * minus_t;
    derivativeShapeFunctions(8, 1) = -0.25 * mid_r * minus_t;
    derivativeShapeFunctions(8, 2) = -0.25 * mid_r * minus_s;

//    shapeFunctions[9] =  0.25 *  plus_r *   mid_s * minus_t;
    derivativeShapeFunctions(9, 0) = +0.25 * mid_s * minus_t;
    derivativeShapeFunctions(9, 1) = -0.5 * s * plus_r * minus_t;
    derivativeShapeFunctions(9, 2) = -0.25 * plus_r * mid_s;

//    shapeFunctions[10]=  0.25 *   mid_r *  plus_s * minus_t;
    derivativeShapeFunctions(10, 0) = -0.5 * r * plus_s * minus_t;
    derivativeShapeFunctions(10, 1) = 0.25 * mid_r * minus_t;
    derivativeShapeFunctions(10, 2) = -0.25 * mid_r * plus_s;

//    shapeFunctions[11]=  0.25 * minus_r *   mid_s * minus_t;
    derivativeShapeFunctions(11, 0) = -0.25 * mid_s * minus_t;
    derivativeShapeFunctions(11, 1) = -0.5 * s * minus_r * minus_t;
    derivativeShapeFunctions(11, 2) = -0.25 * minus_r * mid_s;

//    shapeFunctions[12]=  0.25 * minus_r * minus_s *   mid_t;
    derivativeShapeFunctions(12, 0) = -0.25 * minus_s * mid_t;
    derivativeShapeFunctions(12, 1) = -0.25 * minus_r * mid_t;
    derivativeShapeFunctions(12, 2) = -0.5 * t * minus_r * minus_s;

//    shapeFunctions[13]=  0.25 *  plus_r * minus_s *   mid_t;
    derivativeShapeFunctions(13, 0) = 0.25 * minus_s * mid_t;
    derivativeShapeFunctions(13, 1) = -0.25 * plus_r * mid_t;
    derivativeShapeFunctions(13, 2) = -0.5 * t * plus_r * minus_s;

//    shapeFunctions[14]=  0.25 *  plus_r *  plus_s *   mid_t;
    derivativeShapeFunctions(14, 0) = 0.25 * plus_s * mid_t;
    derivativeShapeFunctions(14, 1) = 0.25 * plus_r * mid_t;
    derivativeShapeFunctions(14, 2) = -0.5 * t * plus_r * plus_s;

//    shapeFunctions[15]=  0.25 * minus_r *  plus_s *   mid_t;
    derivativeShapeFunctions(15, 0) = -0.25 * plus_s * mid_t;
    derivativeShapeFunctions(15, 1) = 0.25 * minus_r * mid_t;
    derivativeShapeFunctions(15, 2) = -0.5 * t * minus_r * plus_s;

//    shapeFunctions[16]=  0.25 *   mid_r * minus_s *  plus_t;
    derivativeShapeFunctions(16, 0) = -0.5 * r * minus_s * plus_t;
    derivativeShapeFunctions(16, 1) = -0.25 * mid_r * plus_t;
    derivativeShapeFunctions(16, 2) = 0.25 * mid_r * minus_s;

//    shapeFunctions[17]=  0.25 *  plus_r *   mid_s *  plus_t;
    derivativeShapeFunctions(17, 0) = 0.25 * mid_s * plus_t;
    derivativeShapeFunctions(17, 1) = -0.5 * s * plus_r * plus_t;
    derivativeShapeFunctions(17, 2) = 0.25 * plus_r * mid_s;

//    shapeFunctions[18]=  0.25 *   mid_r *  plus_s *  plus_t;
    derivativeShapeFunctions(18, 0) = -0.5 * r * plus_s * plus_t;
    derivativeShapeFunctions(18, 1) = 0.25 * mid_r * plus_t;
    derivativeShapeFunctions(18, 2) = 0.25 * mid_r * plus_s;

//    shapeFunctions[19]=  0.25 * minus_r *   mid_s *  plus_t;
    derivativeShapeFunctions(19, 0) = -0.25 * mid_s * plus_t;
    derivativeShapeFunctions(19, 1) = -0.5 * s * minus_r * plus_t;
    derivativeShapeFunctions(19, 2) = 0.25 * minus_r * mid_s;

    return derivativeShapeFunctions;
}

Eigen::Matrix<double, 3, 1> NodeCoordinatesBrickSpectralOrder2(int rNodeIndex)
{
    const int d = 3;
    const int dxd = d * d;

    assert(rNodeIndex >= 0);
    assert(rNodeIndex < d * d * d);

    double cX = ShapeFunctions1D::NodeCoordinatesTrussOrder2(rNodeIndex % d)(0, 0);
    double cY = ShapeFunctions1D::NodeCoordinatesTrussOrder2((rNodeIndex % dxd) / d)(0, 0);
    double cZ = ShapeFunctions1D::NodeCoordinatesTrussOrder2(rNodeIndex / dxd)(0, 0);

    return Eigen::Vector3d(cX, cY, cZ);
}

Eigen::Matrix<double, 27, 1> ShapeFunctionsBrickSpectralOrder2(const Eigen::VectorXd& rCoordinates)
{
    const int d = 3;
    const int dxd = d * d;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, 1, 1>& cZ = rCoordinates.row(2);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussOrder2(cX);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussOrder2(cY);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dz = ShapeFunctions1D::ShapeFunctionsTrussOrder2(cZ);

    Eigen::Matrix<double, d, d> shapeFunctionsMatrix1 = shapeFunctions1Dx * shapeFunctions1Dy.transpose();
    Eigen::Matrix<double, d, dxd> shapeFunctionsMatrix;
    for (int countZ = 0; countZ < d; countZ++)
    {
        shapeFunctionsMatrix.block(0, countZ * d, d, d) = shapeFunctions1Dz(countZ) * shapeFunctionsMatrix1;
    }
    return Eigen::Map<Eigen::Matrix<double, d * dxd, 1>>(shapeFunctionsMatrix.data(), d * dxd);
}

Eigen::Matrix<double, 27, 3> DerivativeShapeFunctionsBrickSpectralOrder2(const Eigen::VectorXd& rCoordinates)
{
    const int d = 3;
    const int dxd = d * d;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, 1, 1>& cZ = rCoordinates.row(2);

    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussOrder2(cX);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussOrder2(cY);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dz = ShapeFunctions1D::ShapeFunctionsTrussOrder2(cZ);

    const Eigen::Matrix<double, d, 1>& derShapeFunctions1Dx = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(cX);
    const Eigen::Matrix<double, d, 1>& derShapeFunctions1Dy = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(cY);
    const Eigen::Matrix<double, d, 1>& derShapeFunctions1Dz = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(cZ);

    Eigen::Matrix<double, dxd * d, 3> derivativeShapeFunctions;

    int theNode(0);
    for (int countzNode = 0; countzNode < d; countzNode++)
        for (int countyNode = 0; countyNode < d; countyNode++)
            for (int countxNode = 0; countxNode < d; countxNode++, theNode++)
            {
                //this can still be optimized, since the calculation of the product of the variables that are not derived is performed several times
                derivativeShapeFunctions(theNode, 0) = derShapeFunctions1Dx(countxNode) * shapeFunctions1Dy(countyNode) * shapeFunctions1Dz(countzNode);
                derivativeShapeFunctions(theNode, 1) = derShapeFunctions1Dy(countyNode) * shapeFunctions1Dz(countzNode) * shapeFunctions1Dx(countxNode);
                derivativeShapeFunctions(theNode, 2) = derShapeFunctions1Dz(countzNode) * shapeFunctions1Dx(countxNode) * shapeFunctions1Dy(countyNode);
            }
    /*   Eigen::Matrix<double,27, 1>  shapeOrig = ShapeFunctionsBrickSpectralOrder2(rCoordinates);
     Eigen::Matrix<double, dxd*d, 3> derivativeShapeFunctionsCDF;
     double delta=1e-8;
     for (int der=0; der<3; der++)
     {
     Eigen::VectorXd coordinates(rCoordinates);
     coordinates(der)+=delta;
     Eigen::Matrix<double,27, 1>  shapeCDF = ShapeFunctionsBrickSpectralOrder2(coordinates);
     derivativeShapeFunctionsCDF.col(der) = 1./delta*(shapeCDF-shapeOrig);
     }
     std::cout << "DerShapeFunctions " << std::endl;
     std::cout << derivativeShapeFunctions << std::endl;
     std::cout << "DerShapeFunctions CDF" << std::endl;
     std::cout << derivativeShapeFunctionsCDF << std::endl;
     std::cout << "diff" << std::endl;
     std::cout << derivativeShapeFunctionsCDF-derivativeShapeFunctions << std::endl;
     */
    return derivativeShapeFunctions;
}

Eigen::Matrix<double, 3, 1> NodeCoordinatesBrickSpectralOrder3(int rNodeIndex)
{
    const int d = 4;
    const int dxd = d * d;

    assert(rNodeIndex >= 0);
    assert(rNodeIndex < d * d * d);

    double cX = ShapeFunctions1D::NodeCoordinatesTrussOrder3(rNodeIndex % d)(0, 0);
    double cY = ShapeFunctions1D::NodeCoordinatesTrussOrder3((rNodeIndex % dxd) / d)(0, 0);
    double cZ = ShapeFunctions1D::NodeCoordinatesTrussOrder3(rNodeIndex / dxd)(0, 0);

    return Eigen::Vector3d(cX, cY, cZ);
}

Eigen::Matrix<double, 64, 1> ShapeFunctionsBrickSpectralOrder3(const Eigen::VectorXd& rCoordinates)
{
    const int d = 4;
    const int dxd = d * d;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, 1, 1>& cZ = rCoordinates.row(2);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussOrder3(cX);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussOrder3(cY);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dz = ShapeFunctions1D::ShapeFunctionsTrussOrder3(cZ);

    Eigen::Matrix<double, d, d> shapeFunctionsMatrix1 = shapeFunctions1Dx * shapeFunctions1Dy.transpose();
    Eigen::Matrix<double, d, dxd> shapeFunctionsMatrix;
    for (int countZ = 0; countZ < d; countZ++)
    {
        shapeFunctionsMatrix.block(0, countZ * d, d, d) = shapeFunctions1Dz(countZ) * shapeFunctionsMatrix1;
    }
    return Eigen::Map<Eigen::Matrix<double, d * dxd, 1>>(shapeFunctionsMatrix.data(), d * dxd);
}

Eigen::Matrix<double, 64, 3> DerivativeShapeFunctionsBrickSpectralOrder3(const Eigen::VectorXd& rCoordinates)
{
    const int d = 4;
    const int dxd = d * d;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, 1, 1>& cZ = rCoordinates.row(2);

    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussOrder3(cX);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussOrder3(cY);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dz = ShapeFunctions1D::ShapeFunctionsTrussOrder3(cZ);

    const Eigen::Matrix<double, d, 1>& derShapeFunctions1Dx = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder3(cX);
    const Eigen::Matrix<double, d, 1>& derShapeFunctions1Dy = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder3(cY);
    const Eigen::Matrix<double, d, 1>& derShapeFunctions1Dz = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder3(cZ);

    Eigen::Matrix<double, dxd * d, 3> derivativeShapeFunctions;

    int theNode(0);
    for (int countzNode = 0; countzNode < d; countzNode++)
        for (int countyNode = 0; countyNode < d; countyNode++)
            for (int countxNode = 0; countxNode < d; countxNode++, theNode++)
            {
                //this can still be optimized, since the calculation of the product of the variables that are not derived is performed several times
                derivativeShapeFunctions(theNode, 0) = derShapeFunctions1Dx(countxNode) * shapeFunctions1Dy(countyNode) * shapeFunctions1Dz(countzNode);
                derivativeShapeFunctions(theNode, 1) = derShapeFunctions1Dy(countyNode) * shapeFunctions1Dz(countzNode) * shapeFunctions1Dx(countxNode);
                derivativeShapeFunctions(theNode, 2) = derShapeFunctions1Dz(countzNode) * shapeFunctions1Dx(countxNode) * shapeFunctions1Dy(countyNode);
            }
    /*   Eigen::Matrix<double,64, 1>  shapeOrig = ShapeFunctionsBrickSpectralOrder3(rCoordinates);
     Eigen::Matrix<double, dxd*d, 3> derivativeShapeFunctionsCDF;
     double delta=1e-8;
     for (int der=0; der<3; der++)
     {
     Eigen::VectorXd coordinates(rCoordinates);
     coordinates(der)+=delta;
     Eigen::Matrix<double,64, 1>  shapeCDF = ShapeFunctionsBrickSpectralOrder3(coordinates);
     derivativeShapeFunctionsCDF.col(der) = 1./delta*(shapeCDF-shapeOrig);
     }
     std::cout << "DerShapeFunctions " << std::endl;
     std::cout << derivativeShapeFunctions << std::endl;
     std::cout << "DerShapeFunctions CDF" << std::endl;
     std::cout << derivativeShapeFunctionsCDF << std::endl;
     std::cout << "diff" << std::endl;
     std::cout << derivativeShapeFunctionsCDF-derivativeShapeFunctions << std::endl;
     */
    return derivativeShapeFunctions;
}

Eigen::Matrix<double, 3, 1> NodeCoordinatesBrickSpectralOrder4(int rNodeIndex)
{
    const int d = 5;
    const int dxd = d * d;

    assert(rNodeIndex >= 0);
    assert(rNodeIndex < d * d * d);

    double cX = ShapeFunctions1D::NodeCoordinatesTrussOrder4(rNodeIndex % d)(0, 0);
    double cY = ShapeFunctions1D::NodeCoordinatesTrussOrder4((rNodeIndex % dxd) / d)(0, 0);
    double cZ = ShapeFunctions1D::NodeCoordinatesTrussOrder4(rNodeIndex / dxd)(0, 0);

    return Eigen::Vector3d(cX, cY, cZ);
}

Eigen::Matrix<double, 125, 1> ShapeFunctionsBrickSpectralOrder4(const Eigen::VectorXd& rCoordinates)
{
    const int d = 5;
    const int dxd = d * d;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, 1, 1>& cZ = rCoordinates.row(2);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussOrder4(cX);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussOrder4(cY);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dz = ShapeFunctions1D::ShapeFunctionsTrussOrder4(cZ);

    Eigen::Matrix<double, d, d> shapeFunctionsMatrix1 = shapeFunctions1Dx * shapeFunctions1Dy.transpose();
    Eigen::Matrix<double, d, dxd> shapeFunctionsMatrix;
    for (int countZ = 0; countZ < d; countZ++)
    {
        shapeFunctionsMatrix.block(0, countZ * d, d, d) = shapeFunctions1Dz(countZ) * shapeFunctionsMatrix1;
    }
    return Eigen::Map<Eigen::Matrix<double, d * dxd, 1>>(shapeFunctionsMatrix.data(), d * dxd);
}

Eigen::Matrix<double, 125, 3> DerivativeShapeFunctionsBrickSpectralOrder4(const Eigen::VectorXd& rCoordinates)
{
    const int d = 5;
    const int dxd = d * d;

    const Eigen::Matrix<double, 1, 1>& cX = rCoordinates.row(0);
    const Eigen::Matrix<double, 1, 1>& cY = rCoordinates.row(1);
    const Eigen::Matrix<double, 1, 1>& cZ = rCoordinates.row(2);

    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dx = ShapeFunctions1D::ShapeFunctionsTrussOrder4(cX);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dy = ShapeFunctions1D::ShapeFunctionsTrussOrder4(cY);
    const Eigen::Matrix<double, d, 1>& shapeFunctions1Dz = ShapeFunctions1D::ShapeFunctionsTrussOrder4(cZ);

    const Eigen::Matrix<double, d, 1>& derShapeFunctions1Dx = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder4(cX);
    const Eigen::Matrix<double, d, 1>& derShapeFunctions1Dy = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder4(cY);
    const Eigen::Matrix<double, d, 1>& derShapeFunctions1Dz = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder4(cZ);

    Eigen::Matrix<double, dxd * d, 3> derivativeShapeFunctions;

    int theNode(0);
    for (int countzNode = 0; countzNode < d; countzNode++)
        for (int countyNode = 0; countyNode < d; countyNode++)
            for (int countxNode = 0; countxNode < d; countxNode++, theNode++)
            {
                //this can still be optimized, since the calculation of the product of the variables that are not derived is performed several times
                derivativeShapeFunctions(theNode, 0) = derShapeFunctions1Dx(countxNode) * shapeFunctions1Dy(countyNode) * shapeFunctions1Dz(countzNode);
                derivativeShapeFunctions(theNode, 1) = derShapeFunctions1Dy(countyNode) * shapeFunctions1Dz(countzNode) * shapeFunctions1Dx(countxNode);
                derivativeShapeFunctions(theNode, 2) = derShapeFunctions1Dz(countzNode) * shapeFunctions1Dx(countxNode) * shapeFunctions1Dy(countyNode);
            }
    /*   Eigen::Matrix<double,125, 1>  shapeOrig = ShapeFunctionsBrickSpectralOrder4(rCoordinates);
     Eigen::Matrix<double, dxd*d, 3> derivativeShapeFunctionsCDF;
     double delta=1e-8;
     for (int der=0; der<3; der++)
     {
     Eigen::VectorXd coordinates(rCoordinates);
     coordinates(der)+=delta;
     Eigen::Matrix<double,125, 1>  shapeCDF = ShapeFunctionsBrickSpectralOrder4(coordinates);
     derivativeShapeFunctionsCDF.col(der) = 1./delta*(shapeCDF-shapeOrig);
     }
     std::cout << "DerShapeFunctions " << std::endl;
     std::cout << derivativeShapeFunctions << std::endl;
     std::cout << "DerShapeFunctions CDF" << std::endl;
     std::cout << derivativeShapeFunctionsCDF << std::endl;
     std::cout << "diff" << std::endl;
     std::cout << derivativeShapeFunctionsCDF-derivativeShapeFunctions << std::endl;
     */
    return derivativeShapeFunctions;
}


Eigen::Matrix<double, 3, 1> NodeCoordinatesPrismOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {

    case 0: return Eigen::Vector3d( 0., 0., -1.);
    case 1: return Eigen::Vector3d( 1., 0., -1.);
    case 2: return Eigen::Vector3d( 0., 1., -1.);
    case 3: return Eigen::Vector3d( 0., 0.,  1.);
    case 4: return Eigen::Vector3d( 1., 0.,  1.);
    case 5: return Eigen::Vector3d( 0., 1.,  1.);
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..5)");
    }
}

Eigen::Matrix<double, 6, 1> ShapeFunctionsPrismOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 6, 1> shapeFunctions;

    double triangle0 = 1. - rCoordinates[0] - rCoordinates[1];
    double triangle1 = rCoordinates[0];
    double triangle2 = rCoordinates[1];

    shapeFunctions[0] = 0.5 * triangle0 * (1. - rCoordinates[2]);
    shapeFunctions[1] = 0.5 * triangle1 * (1. - rCoordinates[2]);
    shapeFunctions[2] = 0.5 * triangle2 * (1. - rCoordinates[2]);
    shapeFunctions[3] = 0.5 * triangle0 * (1. + rCoordinates[2]);
    shapeFunctions[4] = 0.5 * triangle1 * (1. + rCoordinates[2]);
    shapeFunctions[5] = 0.5 * triangle2 * (1. + rCoordinates[2]);
    return shapeFunctions;
}

Eigen::Matrix<double, 6, 3> DerivativeShapeFunctionsPrismOrder1(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 6, 3> derivativeShapeFunctions;

    derivativeShapeFunctions(0, 0) = -0.5 * (1. - rCoordinates[2]);
    derivativeShapeFunctions(0, 1) = -0.5 * (1. - rCoordinates[2]);
    derivativeShapeFunctions(0, 2) = -0.5 * (1. - rCoordinates[0] - rCoordinates[1]);

    derivativeShapeFunctions(1, 0) =  0.5 * (1. - rCoordinates[2]);
    derivativeShapeFunctions(1, 1) =  0.;
    derivativeShapeFunctions(1, 2) = -0.5 * rCoordinates[0];

    derivativeShapeFunctions(2, 0) =  0.;
    derivativeShapeFunctions(2, 1) =  0.5 * (1. - rCoordinates[2]);
    derivativeShapeFunctions(2, 2) = -0.5 * rCoordinates[1];

    derivativeShapeFunctions(3, 0) = -0.5 * (1. + rCoordinates[2]);
    derivativeShapeFunctions(3, 1) = -0.5 * (1. + rCoordinates[2]);
    derivativeShapeFunctions(3, 2) =  0.5 * (1. - rCoordinates[0] - rCoordinates[1]);

    derivativeShapeFunctions(4, 0) =  0.5 * (1. + rCoordinates[2]);
    derivativeShapeFunctions(4, 1) =  0.;
    derivativeShapeFunctions(4, 2) =  0.5 * rCoordinates[0];

    derivativeShapeFunctions(5, 0) =  0.;
    derivativeShapeFunctions(5, 1) =  0.5 * (1. + rCoordinates[2]);
    derivativeShapeFunctions(5, 2) =  0.5 * rCoordinates[1];

    return derivativeShapeFunctions;
}

Eigen::Matrix<double, 3, 1> NodeCoordinatesPrismOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {

    case 0: return Eigen::Vector3d( 0., 0., -1.);
    case 1: return Eigen::Vector3d( 1., 0., -1.);
    case 2: return Eigen::Vector3d( 0., 1., -1.);
    case 3: return Eigen::Vector3d( 0., 0.,  1.);
    case 4: return Eigen::Vector3d( 1., 0.,  1.);
    case 5: return Eigen::Vector3d( 0., 1.,  1.);

    case 6: return Eigen::Vector3d( .5, 0., -1.);
    case 7: return Eigen::Vector3d( 0., .5, -1.);
    case 8: return Eigen::Vector3d( 0., 0.,  0.);
    case 9: return Eigen::Vector3d( .5, .5, -1.);
    case 10: return Eigen::Vector3d( 1., 0.,  0.);
    case 11: return Eigen::Vector3d( 0,  1.,  0.);

    case 12: return Eigen::Vector3d( .5, 0.,  1.);
    case 13: return Eigen::Vector3d( 0., .5,  1.);
    case 14: return Eigen::Vector3d( .5, .5,  1.);

    case 15: return Eigen::Vector3d( .5, 0.,  0.);
    case 16: return Eigen::Vector3d( 0., .5,  0.);
    case 17: return Eigen::Vector3d( .5, .5,  0.);

    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "node index out of range (0..14)");
    }
}

Eigen::Matrix<double, 18, 1> ShapeFunctionsPrismOrder2(const Eigen::VectorXd& rCoordinates)
{
    Eigen::Matrix<double, 18, 1> shapeFunctions;

    auto NTriangle2 = ShapeFunctions2D::ShapeFunctionsTriangleOrder2(Eigen::Vector2d(rCoordinates[0], rCoordinates[1]));
    auto NTruss2 = ShapeFunctions1D::ShapeFunctionsTrussOrder2(Eigen::Matrix<double, 1, 1>::Constant(rCoordinates[2]));

    shapeFunctions[0]  = NTriangle2[0] * NTruss2[0];
    shapeFunctions[1]  = NTriangle2[1] * NTruss2[0];
    shapeFunctions[2]  = NTriangle2[2] * NTruss2[0];
    shapeFunctions[3]  = NTriangle2[0] * NTruss2[2];
    shapeFunctions[4]  = NTriangle2[1] * NTruss2[2];
    shapeFunctions[5]  = NTriangle2[2] * NTruss2[2];
    shapeFunctions[6]  = NTriangle2[3] * NTruss2[0];
    shapeFunctions[7]  = NTriangle2[5] * NTruss2[0];
    shapeFunctions[8]  = NTriangle2[0] * NTruss2[1];
    shapeFunctions[9]  = NTriangle2[4] * NTruss2[0];
    shapeFunctions[10] = NTriangle2[1] * NTruss2[1];
    shapeFunctions[11] = NTriangle2[2] * NTruss2[1];
    shapeFunctions[12] = NTriangle2[3] * NTruss2[2];
    shapeFunctions[13] = NTriangle2[5] * NTruss2[2];
    shapeFunctions[14] = NTriangle2[4] * NTruss2[2];
    shapeFunctions[15] = NTriangle2[3] * NTruss2[1];
    shapeFunctions[16] = NTriangle2[5] * NTruss2[1];
    shapeFunctions[17] = NTriangle2[4] * NTruss2[1];

    return shapeFunctions;
}

Eigen::Matrix<double, 18, 3> DerivativeShapeFunctionsPrismOrder2(const Eigen::VectorXd& rCoordinates)
{

    auto NTriangle2 = ShapeFunctions2D::ShapeFunctionsTriangleOrder2(Eigen::Vector2d(rCoordinates[0], rCoordinates[1]));
    auto dNTriangle2 = ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder2(Eigen::Vector2d(rCoordinates[0], rCoordinates[1]));

    auto NTruss2 = ShapeFunctions1D::ShapeFunctionsTrussOrder2(Eigen::Matrix<double, 1, 1>::Constant(rCoordinates[2]));
    auto dNTruss2 = ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(Eigen::Matrix<double, 1, 1>::Constant(rCoordinates[2]));

    Eigen::Matrix<double, 18, 3> derivativeShapeFunctions;

    derivativeShapeFunctions(0, 0)  = dNTriangle2( 0, 0) *  NTruss2[0];
    derivativeShapeFunctions(0, 1)  = dNTriangle2( 0, 1) *  NTruss2[0];
    derivativeShapeFunctions(0, 2)  =  NTriangle2[0] * dNTruss2[0];
    
    derivativeShapeFunctions(1, 0)  = dNTriangle2( 1, 0) *  NTruss2[0];
    derivativeShapeFunctions(1, 1)  = dNTriangle2( 1, 1) *  NTruss2[0];
    derivativeShapeFunctions(1, 2)  =  NTriangle2[1] * dNTruss2[0];
    
    derivativeShapeFunctions(2, 0)  = dNTriangle2( 2, 0) *  NTruss2[0];
    derivativeShapeFunctions(2, 1)  = dNTriangle2( 2, 1) *  NTruss2[0];
    derivativeShapeFunctions(2, 2)  =  NTriangle2[2] * dNTruss2[0];
    
    derivativeShapeFunctions(3, 0)  = dNTriangle2( 0, 0) *  NTruss2[2];
    derivativeShapeFunctions(3, 1)  = dNTriangle2( 0, 1) *  NTruss2[2];
    derivativeShapeFunctions(3, 2)  =  NTriangle2[0] * dNTruss2[2];
    
    derivativeShapeFunctions(4, 0)  = dNTriangle2( 1, 0) *  NTruss2[2];
    derivativeShapeFunctions(4, 1)  = dNTriangle2( 1, 1) *  NTruss2[2];
    derivativeShapeFunctions(4, 2)  =  NTriangle2[1] * dNTruss2[2];
    
    derivativeShapeFunctions(5, 0)  = dNTriangle2( 2, 0) *  NTruss2[2];
    derivativeShapeFunctions(5, 1)  = dNTriangle2( 2, 1) *  NTruss2[2];
    derivativeShapeFunctions(5, 2)  =  NTriangle2[2] * dNTruss2[2];
    
    derivativeShapeFunctions(6, 0)  = dNTriangle2( 3, 0) *  NTruss2[0];
    derivativeShapeFunctions(6, 1)  = dNTriangle2( 3, 1) *  NTruss2[0];
    derivativeShapeFunctions(6, 2)  =  NTriangle2[3] * dNTruss2[0];
    
    derivativeShapeFunctions(7, 0)  = dNTriangle2( 5, 0) *  NTruss2[0];
    derivativeShapeFunctions(7, 1)  = dNTriangle2( 5, 1) *  NTruss2[0];
    derivativeShapeFunctions(7, 2)  =  NTriangle2[5] * dNTruss2[0];
    
    derivativeShapeFunctions(8, 0)  = dNTriangle2( 0, 0) *  NTruss2[1];
    derivativeShapeFunctions(8, 1)  = dNTriangle2( 0, 1) *  NTruss2[1];
    derivativeShapeFunctions(8, 2)  =  NTriangle2[0] * dNTruss2[1];
    
    derivativeShapeFunctions(9, 0)  = dNTriangle2( 4, 0) *  NTruss2[0];
    derivativeShapeFunctions(9, 1)  = dNTriangle2( 4, 1) *  NTruss2[0];
    derivativeShapeFunctions(9, 2)  =  NTriangle2[4] * dNTruss2[0];
    
    derivativeShapeFunctions(10, 0) = dNTriangle2( 1, 0) *  NTruss2[1];
    derivativeShapeFunctions(10, 1) = dNTriangle2( 1, 1) *  NTruss2[1];
    derivativeShapeFunctions(10, 2) =  NTriangle2[1] * dNTruss2[1];
    
    derivativeShapeFunctions(11, 0) = dNTriangle2( 2, 0) *  NTruss2[1];
    derivativeShapeFunctions(11, 1) = dNTriangle2( 2, 1) *  NTruss2[1];
    derivativeShapeFunctions(11, 2) =  NTriangle2[2] * dNTruss2[1];
    
    derivativeShapeFunctions(12, 0) = dNTriangle2( 3, 0) *  NTruss2[2];
    derivativeShapeFunctions(12, 1) = dNTriangle2( 3, 1) *  NTruss2[2];
    derivativeShapeFunctions(12, 2) =  NTriangle2[3] * dNTruss2[2];
    
    derivativeShapeFunctions(13, 0) = dNTriangle2( 5, 0) *  NTruss2[2];
    derivativeShapeFunctions(13, 1) = dNTriangle2( 5, 1) *  NTruss2[2];
    derivativeShapeFunctions(13, 2) =  NTriangle2[5] * dNTruss2[2];
    
    derivativeShapeFunctions(14, 0) = dNTriangle2( 4, 0) *  NTruss2[2];
    derivativeShapeFunctions(14, 1) = dNTriangle2( 4, 1) *  NTruss2[2];
    derivativeShapeFunctions(14, 2) =  NTriangle2[4] * dNTruss2[2];

    derivativeShapeFunctions(15, 0) = dNTriangle2( 3, 0) *  NTruss2[1];
    derivativeShapeFunctions(15, 1) = dNTriangle2( 3, 1) *  NTruss2[1];
    derivativeShapeFunctions(15, 2) =  NTriangle2[3] * dNTruss2[1];

    derivativeShapeFunctions(16, 0) = dNTriangle2( 5, 0) *  NTruss2[1];
    derivativeShapeFunctions(16, 1) = dNTriangle2( 5, 1) *  NTruss2[1];
    derivativeShapeFunctions(16, 2) =  NTriangle2[5] * dNTruss2[1];

    derivativeShapeFunctions(17, 0) = dNTriangle2( 4, 0) *  NTruss2[1];
    derivativeShapeFunctions(17, 1) = dNTriangle2( 4, 1) *  NTruss2[1];
    derivativeShapeFunctions(17, 2) =  NTriangle2[4] * dNTruss2[1];
    
    return derivativeShapeFunctions;
}

}


namespace ShapeFunctionsInterface2D // interval -1 to 1
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::MatrixXd NodeCoordinatesInterface2dOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Vector2d(-1, -1);
    case 1: return Eigen::Vector2d(+1, -1);
    case 2: return Eigen::Vector2d(+1, +1);
    case 3: return Eigen::Vector2d(-1, +1);
    default:
        throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t node index out of range (0..3)");
    }
}

Eigen::MatrixXd ShapeFunctionsInterface2dOrder1(const Eigen::VectorXd& rCoordinates)
{
    const double N00 = 0.5 * (1. - rCoordinates[0]);
    const double N01 = 0.5 * (1. + rCoordinates[0]);

    return Eigen::Vector4d(-N00, -N01, N01, N00);
}

Eigen::MatrixXd DerivativeShapeFunctionsInterface2dOrder1(const Eigen::VectorXd& rCoordinates)
{
    // this interface element does not need any shape function derivatives
    return Eigen::Matrix2d::Zero();
}

Eigen::MatrixXd NodeCoordinatesInterface2dOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Vector2d(-1., -1.);
    case 1: return Eigen::Vector2d(+0., -1);
    case 2: return Eigen::Vector2d(+1, -1);
    case 3: return Eigen::Vector2d(+1, +1);
    case 4: return Eigen::Vector2d(0, +1);
    case 5: return Eigen::Vector2d(-1, +1);
    default:
        throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t node index out of range (0..5)");
    }
}

Eigen::MatrixXd ShapeFunctionsInterface2dOrder2(const Eigen::VectorXd& rCoordinates)
{
    const double xi = rCoordinates(0,0);

    const double N00 = 0.5 * xi * (xi - 1.);
    const double N01 = (1. - xi * xi);
    const double N02 = 0.5 * xi * (xi + 1.);

    return (Eigen::MatrixXd(6,1) << N00, N01, N02, -N02, -N01, -N00).finished();
}

Eigen::MatrixXd DerivativeShapeFunctionsInterface2dOrder2(const Eigen::VectorXd& rCoordinates)
{
    // this interface element does not need any shape function derivatives
    return Eigen::Matrix2d::Zero();
}


}

namespace ShapeFunctionsInterface3D // interval -1 to 1
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::MatrixXd NodeCoordinatesInterface3dOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0: return Eigen::Vector3d(-1, -1, 0);
    case 1: return Eigen::Vector3d(+1, -1, 0);
    case 2: return Eigen::Vector3d(+1, +1, 0);
    case 3: return Eigen::Vector3d(-1, +1, 0);
    default:
        throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t node index out of range (0..3)");
    }
}

Eigen::MatrixXd ShapeFunctionsInterface3dOrder1(const Eigen::VectorXd& rCoordinates)
{
    const double N00 = 0.5 * (1. - rCoordinates[0]);
    const double N01 = 0.5 * (1. + rCoordinates[0]);

    return Eigen::Vector4d(-N00, -N01, N01, N00);
}

Eigen::MatrixXd DerivativeShapeFunctionsInterface3dOrder1(const Eigen::VectorXd& rCoordinates)
{
    // this interface element does not need any shape function derivatives
    return Eigen::Matrix3d::Zero();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}// namespace ShapeFunctionsInterface3D


// In order to maintain equal shape functions on each element the Bézier extraction together with Bernstein polynomials is used (see. Borden et. al. 2011)
namespace ShapeFunctionsIGA
{

/////////////////////////////// BSPLINE ////////////////////////////////////////////////////////////////////////

int FindSpan(double rParameter, int rDegree, const Eigen::VectorXd &rKnots)
{
    int size = rKnots.rows();
    assert(rParameter >= rKnots(0) && rParameter <= rKnots(size - 1));
    int numBasisFuns = size - rDegree - 1;
    if(rParameter == rKnots[numBasisFuns]) return numBasisFuns-1;

    int low  = rDegree;
    int high = numBasisFuns;
    int mid  = (low + high)/2;

    while(rParameter < rKnots[mid] || rParameter >= rKnots[mid+1])
    {
        if(rParameter < rKnots[mid]) high = mid;
        else                         low = mid;

        mid = (low + high)/2;
    }
    return mid;
}

Eigen::VectorXd BasisFunctions(double rParameter, int spanIdx, int rDegree, const Eigen::VectorXd &rKnots)
{
    Eigen::VectorXd rBasisFunctions(rDegree+1);

    rBasisFunctions[0] = 1.;

    Eigen::VectorXd left(rDegree+1);
    Eigen::VectorXd right(rDegree+1);

    for (int j = 1; j <= rDegree; j++)
    {
        left[j]  = rParameter - rKnots[spanIdx + 1 - j];
        right[j] = rKnots[spanIdx + j] - rParameter;
        double saved = 0.;
        for(int r = 0; r < j; r++)
        {
            double temp = rBasisFunctions[r]/(right[r+1] + left[j-r]);
            rBasisFunctions[r] = saved + right[r+1]*temp;
            saved = left[j-r]*temp;
        }
        rBasisFunctions[j] = saved;
    }

    return rBasisFunctions;
}

Eigen::MatrixXd BasisFunctionsAndDerivatives( int der, double rParameter, int spanIdx, int rDegree, const Eigen::VectorXd &rKnots)
{
    Eigen::MatrixXd ndu(rDegree+1, rDegree+1);

    ndu(0,0) = 1.0;

    Eigen::VectorXd left(rDegree + 1);
    Eigen::VectorXd right(rDegree + 1);

    for (int j = 1; j <= rDegree; j++)
    {
        left[j]  = rParameter - rKnots[spanIdx + 1 - j];
        right[j] = rKnots[spanIdx + j] - rParameter;
        double saved = 0.;
        for(int r = 0; r < j; r++)
        {
            ndu(j,r) = right[r+1] + left[j-r]; // lower triange (knot differences)
            double temp = ndu(r, j-1)/ndu(j,r);
            ndu(r,j) = saved + right[r+1]*temp; // upper triangle
            saved = left[j-r]*temp;
        }
        ndu(j,j) = saved;
    }

    Eigen::MatrixXd ders(der+1, rDegree+1);

    for(int j = 0; j <= rDegree; j++) ders(0,j) = ndu(j,rDegree);

    Eigen::MatrixXd a(2, rDegree+1);
    for(int r = 0; r<=rDegree; r++)
    {
        int s1 = 0;
        int s2 = 1;

        a(0,0) = 1.0;

        int j = 0;
        for(int k = 1; k<=der; k++)
        {
            double d = 0.;
            int rk = r-k;
            int pk = rDegree-k;
            if(r >= k)
            {
                a(s2,0) = a(s1,0)/ndu(pk+1, rk);
                d = a(s2,0)*ndu(rk,pk);
            }
            int j1 = 0;
            if(rk >= -1) j1 = 1;
            else         j1 = -rk;

            int j2 = 0;
            if(r-1 <= pk) j2 = k-1;
            else          j2 = rDegree-r;

            for(j = j1; j <= j2; j++)
            {
                a(s2, j) = (a(s1,j) - a(s1, j-1))/ndu(pk+1, rk+j);
                d += a(s2,j)*ndu(rk+j, pk);
            }

            if(r <= pk)
            {
                a(s2,k) = -a(s1,k-1)/ndu(pk+1,r);
                d += a(s2,k)*ndu(r,pk);
            }

            ders(k,r) = d;
            j = s1;
            s1 = s2;
            s2 = j;
        }
    }

    int r = rDegree;
    for(int k = 1; k <= der; k++)
    {
        for (int j = 0; j <= rDegree; j++) ders(k,j) *=r;
        r*= (rDegree-k);
    }


    return ders;
}

Eigen::VectorXd BasisFunctionsRat(double rParameter, int spanIdx, int rDegree, const Eigen::VectorXd &rKnots, const Eigen::VectorXd &rWeights)
{
    Eigen::VectorXd rBasisFunctions = BasisFunctions(rParameter, spanIdx, rDegree, rKnots);

    // NURBS specific ...
    double sum = 0.;
    for(int i = 0; i <= rDegree; i++) sum += rBasisFunctions(i)*rWeights(spanIdx - rDegree + i);
    for(int i = 0; i <= rDegree; i++) rBasisFunctions(i) *= (rWeights(spanIdx - rDegree + i)/sum);

    return rBasisFunctions;
}

Eigen::VectorXd BasisFunctionsAndDerivativesRat(int der, double rParameter, int spanIdx, int rDegree, const Eigen::VectorXd &rKnots, const Eigen::VectorXd &rWeights)
{
    if(der < 0 || der > 2)
        throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t der greater than 2 not implemented, possible values 0,1,2!");

    Eigen::MatrixXd ders = BasisFunctionsAndDerivatives(der, rParameter, spanIdx, rDegree, rKnots);

    // NURBS specific ...
    Eigen::VectorXd sum(der+1);
    sum.setZero(der+1);
    for(int numDer = 0; numDer < der; numDer++)
    {
        for(int i = 0; i <= rDegree; i++)
        {
            sum(numDer) += ders(numDer,i)*rWeights(spanIdx - rDegree + i);
        }
    }


    Eigen::VectorXd dersRat(rDegree+1);
    dersRat.setZero(rDegree+1);

    for(int i = 0; i <= rDegree; i++)
    {
        double weight = rWeights(spanIdx - rDegree + i);
        if     (der == 0)
        {
            dersRat(i) = ders(der, i)*weight/sum(0);
        }
        else if(der == 1)
        {
            dersRat(i) = (ders(der, i)*sum(0) - ders(0,i)*sum(1))* weight/(sum(0)*sum(0));
        }
        else
        {
            double sum2 = sum(0)*sum(0);
            dersRat(i) = weight*(ders(der, i)/sum(0) - 2*ders(1,i)*sum(1)/(sum2) - ders(0,i)*sum(2)/(sum2) + 2*ders(0,i)*sum(1)*sum(1)/(sum2*sum(0)));
        }
    }

    return dersRat;
}

Eigen::VectorXd BasisFunctions2DRat(const Eigen::VectorXd &rCoordinates,
                                    const Eigen::Vector2i &rSpanIdx,
                                    const Eigen::Vector2i &rDegree,
                                    const Eigen::VectorXd &rKnotsX,
                                    const Eigen::VectorXd &rKnotsY,
                                    const Eigen::MatrixXd &rWeights)
{
    Eigen::VectorXd xBasis = BasisFunctions(rCoordinates(0), rSpanIdx(0), rDegree(0), rKnotsX);
    Eigen::VectorXd yBasis = BasisFunctions(rCoordinates(1), rSpanIdx(1), rDegree(1), rKnotsY);

    Eigen::VectorXd basis((rDegree(0)+1)*(rDegree(1)+1));
    basis.setZero((rDegree(0)+1)*(rDegree(1)+1));

    double sum = 0.;
    for(int i = 0; i <= rDegree(1); i++)
    {
        for(int j = 0; j <= rDegree(0); j++)
        {
            double weight = rWeights(rSpanIdx(1) - rDegree(1) + i, rSpanIdx(0) - rDegree(0) + j);
            sum += xBasis(j)*yBasis(i)*weight;
        }
    }

    int count = 0;
    for(int i = 0; i <= rDegree(1); i++)
    {
        for(int j = 0; j <= rDegree(0); j++)
        {
            double weight = rWeights(rSpanIdx(1) - rDegree(1) + i, rSpanIdx(0) - rDegree(0) + j);
            basis(count) = xBasis(j)*yBasis(i)*weight/sum;
            count++;
        }
    }

    return basis;
}

Eigen::MatrixXd BasisFunctionsAndDerivatives2DRat(int                    der,
                                                  const Eigen::VectorXd &rCoordinates,
                                                  const Eigen::Vector2i &rSpanIdx,
                                                  const Eigen::Vector2i &rDegree,
                                                  const Eigen::VectorXd &rKnotsX,
                                                  const Eigen::VectorXd &rKnotsY,
                                                  const Eigen::MatrixXd &rWeights)
{
    if(der < 0 || der >2)
        throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t 'der' greater than 2 not implemented, possible values 0,1,2!");

    Eigen::MatrixXd xBasisDer = BasisFunctionsAndDerivatives(der, rCoordinates(0), rSpanIdx(0), rDegree(0), rKnotsX);
    Eigen::MatrixXd yBasisDer = BasisFunctionsAndDerivatives(der, rCoordinates(1), rSpanIdx(1), rDegree(1), rKnotsY);

    Eigen::MatrixXd ders((rDegree(0)+1)*(rDegree(1)+1), der+1);
    ders.setZero((rDegree(0)+1)*(rDegree(1)+1), der+1);


    int num = ( (der+1)*(der+1) + (der+1) )/2;
    Eigen::VectorXd sum(num);
    sum.setZero(num);

    for(int i = 0; i <= rDegree(1); i++)
    {
        for(int j = 0; j <= rDegree(0); j++)
        {
            double weight = rWeights(rSpanIdx(1) - rDegree(1) + i, rSpanIdx(0) - rDegree(0) + j);

            if     (der == 0)
            {
                sum(0) += xBasisDer(0, j)*yBasisDer(0, i)*weight;
            }
            else if(der == 1)
            {
                sum(0) += xBasisDer(0, j)*yBasisDer(0, i)*weight;
                sum(1) += xBasisDer(1, j)*yBasisDer(0, i)*weight;
                sum(2) += xBasisDer(0, j)*yBasisDer(1, i)*weight;
            }
            else
            {
                sum(0) += xBasisDer(0, j)*yBasisDer(0, i)*weight;
                sum(1) += xBasisDer(1, j)*yBasisDer(0, i)*weight;
                sum(2) += xBasisDer(2, j)*yBasisDer(0, i)*weight;

                sum(3) += xBasisDer(0, j)*yBasisDer(1, i)*weight;
                sum(4) += xBasisDer(0, j)*yBasisDer(2, i)*weight;

                sum(5) += xBasisDer(1, j)*yBasisDer(1, i)*weight;
            }
        }
    }

    int count = 0;
    for(int i = 0; i <= rDegree(1); i++)
    {
        for(int j = 0; j <= rDegree(0); j++)
        {
            double weight = rWeights(rSpanIdx(1) - rDegree(1) + i, rSpanIdx(0) - rDegree(0) + j);

            if     (der == 0)
            {
                ders(count,0) = xBasisDer(0,j)*yBasisDer(0,i)*weight/sum(0);
                count++;
            }
            else if(der == 1)
            {
                ders(count,0) = (xBasisDer(1,j)*yBasisDer(0,i)*sum(0)
                               - xBasisDer(0,j)*yBasisDer(0,i)*sum(1)) * weight/(sum(0)*sum(0));
                ders(count,1) = (xBasisDer(0,j)*yBasisDer(1,i)*sum(0)
                               - xBasisDer(0,j)*yBasisDer(0,i)*sum(2)) * weight/(sum(0)*sum(0));
                count++;
            }
            else
            {
                ders(count,0) =(xBasisDer(2,j)*yBasisDer(0,i)/sum(0)
                            - 2*xBasisDer(1,j)*yBasisDer(0,i)*sum(1)/(sum(0)*sum(0))
                              - xBasisDer(0,j)*yBasisDer(0,i)*sum(2)/(sum(0)*sum(0))
                            + 2*xBasisDer(0,j)*yBasisDer(0,i)*sum(1)*sum(1)/(sum(0)*sum(0)*sum(0)) ) * weight;

                ders(count,1) = (xBasisDer(0,j)*yBasisDer(2,i)/sum(0)
                             - 2*xBasisDer(0,j)*yBasisDer(1,i)*sum(3)/(sum(0)*sum(0))
                               - xBasisDer(0,j)*yBasisDer(0,i)*sum(4)/(sum(0)*sum(0))
                             + 2*xBasisDer(0,j)*xBasisDer(0,i)*sum(3)*sum(3)/(sum(0)*sum(0)*sum(0)) ) * weight;

                ders(count,2) = (xBasisDer(1,j)*yBasisDer(1,i)/sum(0)
                               - xBasisDer(1,j)*yBasisDer(0,i)*sum(3)/(sum(0)*sum(0))
                               - xBasisDer(0,j)*yBasisDer(1,i)*sum(1)/(sum(0)*sum(0))
                               - xBasisDer(1,j)*yBasisDer(0,i)*sum(5)*sum(5)/(sum(0)*sum(0))
                             + 2*xBasisDer(0,j)*yBasisDer(0,i)*sum(1)*sum(3)/(sum(0)*sum(0)*sum(0)) ) * weight;
                count++;
            }
        }
    }

    return ders;
}

///////////////////////////////////// BERNSTEIN ////////////////////////////////////////////////////////////////

Eigen::VectorXd Bernstein1DOrder1(double rParameter)
{
    return Bernstein1D(rParameter, 1);
}

Eigen::VectorXd Bernstein1DOrder2(double rParameter)
{
    return Bernstein1D(rParameter, 2);
}

Eigen::VectorXd Bernstein1DOrder3(double rParameter)
{
    return Bernstein1D(rParameter, 3);
}

Eigen::VectorXd Bernstein1DOrder4(double rParameter)
{
    return Bernstein1D(rParameter, 4);
}

Eigen::VectorXd Bernstein1D(double rParameter, int rOrder)
{
    Eigen::VectorXd values(rOrder);

    double paramUnivariate = (rParameter + 1.)/2;
    double u1 = 1. - paramUnivariate;

    for(int j = 1; j < rOrder; j++)
    {
        double saved = 0.;
        for(int k = 0; k < j; k++)
        {
            double temp = values(k);
            values(k) = saved + u1*temp;
            saved = paramUnivariate*temp;
        }
        values(j) = saved;
    }
    return values;
}


Eigen::VectorXd DerivativeBernstein1DOrder1(double rParameter)
{
    return DerivativeBernstein1D(rParameter, 1);
}

Eigen::VectorXd DerivativeBernstein1DOrder2(double rParameter)
{
    return DerivativeBernstein1D(rParameter, 2);
}

Eigen::VectorXd DerivativeBernstein1DOrder3(double rParameter)
{
    return DerivativeBernstein1D(rParameter, 3);
}

Eigen::VectorXd DerivativeBernstein1DOrder4(double rParameter)
{
    return DerivativeBernstein1D(rParameter, 4);
}

Eigen::VectorXd DerivativeBernstein1D(double rParameter, int rOrder)
{
    Eigen::VectorXd values(rOrder);
    double u1 = rOrder-1;

    for(int j = 1; j < rOrder; j++)
    {
        double saved = 0.;
        for(int k = 0; k < j; k++)
        {
            double temp = values(k);
            values(k) = saved + u1*temp;
            saved = u1*temp;
        }
        values(j) = saved;
    }
    return values;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}// namespace ShapeFunctionsIGA1D



}

