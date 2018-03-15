#pragma once

#include <iosfwd>
#include <string>
#include <Eigen/Core>

namespace NuTo
{

class GmshWriter
{
public:
    struct Box2D
    {
        Eigen::Vector2d mStart;
        Eigen::Vector2d mEnd;
    };

    struct Box3D
    {
        Eigen::Vector3d mStart;
        Eigen::Vector3d mEnd;
    };

    struct Cylinder
    {
        double mRadius;
        double mHeight;
    };

    struct Options
    {
        static Options Default(double meshSizeAggregates, double meshSizeContainer, double interfaceThickness = 0.)
        {
            Options opt;
            opt.interfaceThickness = interfaceThickness;
            opt.meshSizeAggregates = meshSizeAggregates;
            opt.meshSizeContainer = meshSizeContainer;
            return opt;
        }

        double interfaceThickness = 0;
        double meshSizeAggregates = 1;
        double meshSizeContainer = 1;

        int meshAlg = 6;
        int meshAlg3D = 4;
        int meshRecombinationAlgorithm = 0;
        int meshOptimize = 2;
        int meshSmoothing = 2;

        std::string physicalGroupMatrix = "Matrix";
        std::string physicalGroupAggregates = "Aggregates";
        std::string physicalGroupInterfaces = "Interfaces";
        std::string additionalOptions = "";
    };

    static void Write(std::ostream& out, Box2D box, const Eigen::MatrixX3d& aggregates, Options opt);
    static void Write(std::ostream& out, Box3D box, const Eigen::MatrixX4d& aggregates, Options opt);
    static void Write(std::ostream& out, Cylinder cylinder, const Eigen::MatrixX4d& aggregates, Options opt);

    static void Write(std::string filename, Box2D box, const Eigen::MatrixX3d& aggregates, Options opt);
    static void Write(std::string filename, Box3D box, const Eigen::MatrixX4d& aggregates, Options opt);
    static void Write(std::string filename, Cylinder cylinder, const Eigen::MatrixX4d& aggregates, Options opt);
};
} /* NuTo */
