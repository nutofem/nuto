#pragma once

#include <iosfwd>
#include <string>
#include <Eigen/Core>

namespace NuTo
{

//! Writes a mesoscale concrete geometry, consisting of a matrix, several aggregates and interfaces around them. The
//! interfaces are created by extruding the aggregate mesh. This special feature is not available in the gmsh3 open
//! cascade kernel. Thus, it requires very long and error prone gmsh code and, thus, IMO, worth wrapping in this class.
//!
//! This class is somehow connected to the concrete geometry. However, for decoupling and easy external use, it will not
//! rely on the other classes in geometryConcrete. So it will not work on a `NuTo::Specimen` and a
//! `NuTo::ParticleHandler`. It will rather define own simplified boxes and take the aggregates as Eigen::Matrix.
//!
//! The 2D sphere format is
//!   x0  y0  r0
//!   x1  y1  r1
//!   x2  y2  r2
//!  ... ... ...
//!
//! The 3D sphere format is
//!   x0  y0  z0  r0
//!   x1  y1  z1  r1
//!   x2  y2  z2  r2
//!  ... ... ... ...
//!
class GmshWriter
{
public:
    //! Simple 2D box geometry
    struct Box2D
    {
        Box2D(double lx, double ly)
            : mStart(Eigen::Vector2d::Zero())
            , mEnd(Eigen::Vector2d(lx, ly))
        {
        }
        Box2D(Eigen::Vector2d start, Eigen::Vector2d end)
            : mStart(start)
            , mEnd(end)
        {
        }

        Eigen::Vector2d mStart;
        Eigen::Vector2d mEnd;
    };

    //! Simple 3D box geometry
    struct Box3D
    {
        Box3D(double lx, double ly, double lz)
            : mStart(Eigen::Vector3d::Zero())
            , mEnd(Eigen::Vector3d(lx, ly, lz))
        {
        }
        Box3D(Eigen::Vector3d start, Eigen::Vector3d end)
            : mStart(start)
            , mEnd(end)
        {
        }

        Eigen::Vector3d mStart;
        Eigen::Vector3d mEnd;
    };

    //! Simple 3D cylinder geometry
    //! @remark The center will be at (0,0,0). You can make that changable. YAGNI otherwise.
    struct Cylinder
    {
        double mRadius;
        double mHeight;
    };

    //! Gmsh options that are defaulted to well working parameters. You can override them here and add arbitrary options
    //! via the member `additionalOptions`. One idea would be to specifiy `additionalOptions =
    //! "Mesh.SecondOrderIncomplete = 1\n;"` to get 8 node quads in a 2D interface, instead of 9 node quads
    struct Options
    {
        Options(double meshSizeAggregates, double meshSizeMatrix, double interfaceThickness = 0.)
            : interfaceThickness(interfaceThickness)
            , meshSizeAggregates(meshSizeAggregates)
            , meshSizeMatrix(meshSizeMatrix)
        {
        }

        double interfaceThickness = 0;
        double meshSizeAggregates = 1;
        double meshSizeMatrix = 1;

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

    // Write to output stream

    static void Write(std::ostream& out, Box2D box, const Eigen::MatrixX3d& aggregates, Options opt);
    static void Write(std::ostream& out, Box3D box, const Eigen::MatrixX4d& aggregates, Options opt);
    static void Write(std::ostream& out, Cylinder cylinder, const Eigen::MatrixX4d& aggregates, Options opt);

    // Write to file

    static void Write(std::string filename, Box2D box, const Eigen::MatrixX3d& aggregates, Options opt);
    static void Write(std::string filename, Box3D box, const Eigen::MatrixX4d& aggregates, Options opt);
    static void Write(std::string filename, Cylinder cylinder, const Eigen::MatrixX4d& aggregates, Options opt);
};
} /* NuTo */
