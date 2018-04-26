#include "nuto/base/Exception.h"
#include "InterpolationBrickQuadratic.h"

using namespace NuTo;

Eigen::Matrix<double, 3, 1> InterpolationBrickQuadratic::NodeCoordinatesBrickOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    // vertices
    case 0:
        return Eigen::Vector3d(-1., -1., -1.);
    case 1:
        return Eigen::Vector3d(1., -1., -1.);
    case 2:
        return Eigen::Vector3d(1., 1., -1.);
    case 3:
        return Eigen::Vector3d(-1., 1., -1.);
    case 4:
        return Eigen::Vector3d(-1., -1., 1.);
    case 5:
        return Eigen::Vector3d(1., -1., 1.);
    case 6:
        return Eigen::Vector3d(1., 1., 1.);
    case 7:
        return Eigen::Vector3d(-1., 1., 1.);

    // "rear" plane
    case 8:
        return Eigen::Vector3d(0., -1., -1.);
    case 9:
        return Eigen::Vector3d(1., 0., -1.);
    case 10:
        return Eigen::Vector3d(0., 1., -1.);
    case 11:
        return Eigen::Vector3d(-1., 0., -1.);

    // "middle" plane
    case 12:
        return Eigen::Vector3d(-1., -1., 0.);
    case 13:
        return Eigen::Vector3d(1., -1., 0.);
    case 14:
        return Eigen::Vector3d(1., 1., 0.);
    case 15:
        return Eigen::Vector3d(-1., 1., 0.);

    // "front" plane
    case 16:
        return Eigen::Vector3d(0., -1., 1.);
    case 17:
        return Eigen::Vector3d(1., 0., 1.);
    case 18:
        return Eigen::Vector3d(0., 1., 1.);
    case 19:
        return Eigen::Vector3d(-1., 0., 1.);

    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "node index out of range (0..19)");
    }
}

Eigen::Matrix<double, 20, 1> InterpolationBrickQuadratic::ShapeFunctionsBrickOrder2(const Eigen::VectorXd& rCoordinates)
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

Eigen::Matrix<double, 20, 3> InterpolationBrickQuadratic::DerivativeShapeFunctionsBrickOrder2(const Eigen::VectorXd& rCoordinates)
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
