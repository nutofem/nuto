#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Core>

using ExactStress = Eigen::Vector3d;
constexpr double load = 10;

namespace NuTo
{
namespace Test
{

class PlateWithHoleAnalytical
{
public:
    // @brief returns the analytical solution (sxx, syy, sxy) taken from
    // https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Plate_with_hole_in_tension
    static ExactStress AnalyticStress(Eigen::VectorXd cartesianCoordinate)
    {
        constexpr double a = 1.0;
        double r = cartesianCoordinate.norm();
        double theta = std::atan(cartesianCoordinate[1] / cartesianCoordinate[0]);

        double cos2t = std::cos(2 * theta);
        double cos4t = std::cos(4 * theta);
        double sin2t = std::sin(2 * theta);
        double sin4t = std::sin(4 * theta);

        double fac1 = (a * a) / (r * r);
        double fac2 = 1.5 * fac1 * fac1;
        ExactStress stress;
        stress[0] = 1. - fac1 * (1.5 * cos2t + cos4t) + fac2 * cos4t;
        stress[1] = -fac1 * (0.5 * cos2t - cos4t) - fac2 * cos4t;
        stress[2] = -fac1 * (0.5 * sin2t + sin4t) + fac2 * sin4t;

        return stress * load;
    }

    static Eigen::VectorXd AnalyticDisplacement(Eigen::Vector2d cartesianCoordinate, double E, double nu)
    {
        double a = 1.0;
        double T = -load;
        double r = cartesianCoordinate.norm();
        double theta = std::atan(cartesianCoordinate[1] / cartesianCoordinate[0]);

        double ct = std::cos(theta);
        double c3t = std::cos(3 * theta);
        double st = std::sin(theta);
        double s3t = std::sin(3 * theta);

        double mu = E / (2. + 2. * nu);
        double k = (3. - nu) / (1. + nu);

        double ux = ((T * a) / (8. * mu)) * ((r / a) * (k + 1.) * ct + ((2. * a) / r) * ((1. + k) * ct + c3t) -
                                             ((2. * a * a * a) / (r * r * r)) * c3t);
        double uy = ((T * a) / (8. * mu)) * ((r / a) * (k - 3.) * st + ((2. * a) / r) * ((1. - k) * st + s3t) -
                                             ((2. * a * a * a) / (r * r * r)) * s3t);

        return Eigen::Vector2d(ux, uy);
    }

    static Eigen::Vector2d Pressure(Eigen::Vector2d cartesianCoordinate, Eigen::Vector2d rN)
    {
        ExactStress s = AnalyticStress(cartesianCoordinate);
        Eigen::Matrix2d stress;
        stress << s[0], s[2], s[2], s[1];
        Eigen::Vector2d pressure = stress * rN;
        std::cout << "Apply pressure " << pressure.transpose() << " at cartesian coordinate "
                  << cartesianCoordinate.transpose() << std::endl;
        return pressure;
    }

    static Eigen::Vector2d PressureRight(Eigen::Vector2d cartesianCoordinate)
    {
        return Pressure(cartesianCoordinate, Eigen::Vector2d::UnitX());
    }

    static Eigen::Vector2d PressureTop(Eigen::Vector2d cartesianCoordinate)
    {
        return Pressure(cartesianCoordinate, Eigen::Vector2d::UnitY());
    }

    static Eigen::Vector2d PressureLeft(Eigen::Vector2d cartesianCoordinate)
    {
        return Pressure(cartesianCoordinate, -Eigen::Vector2d::UnitX());
    }

    static Eigen::Vector2d PressureBottom(Eigen::Vector2d cartesianCoordinate)
    {
        return Pressure(cartesianCoordinate, -Eigen::Vector2d::UnitY());
    }
};
} // namespace Test
} // namespace NuTo
