/*
 * GeometryConcreteAPI.cpp
 *
 *  Created on: 4 Sep 2015
 *      Author: ttitsche
 */
#include "geometryConcrete/GeometryConcrete.h"
#include <boost/filesystem.hpp>
#include <iostream>


void MaxDistanceMesh2D(std::string rGmshFile, double rLX, double rLY, int rNumClasses)
{
    // ***************************************************************************
    // ***********   define the geometry      ************************************

    NuTo::GeometryConcrete geometry;
    geometry.SetSeed(1337);
    geometry.SetSpecimenBox(0, rLX, 0, rLY, 0, rLY);
    geometry.SetGradingCurve(NuTo::GeometryConcrete::B16, rNumClasses);
    geometry.SetParticleVolumeFraction(0.4);
    geometry.SetAbsoluteGrowthRate(0.1);
    geometry.SetInitialTimeBarrier(10);

    geometry.MaximizeParticleDistance(0.75);

    geometry.ExportGmshGeo2D(rGmshFile, 0.75, rLY/2., 0.75);
}

void MaxVolumeFraction3D(std::string rGmshFile, double rLX, double rLY, double rLZ, int rNumClasses)
{
    // ***************************************************************************
    // ***********   define the geometry      ************************************

    NuTo::GeometryConcrete geometry;
    geometry.SetSeed(1337);
    geometry.SetSpecimenBox(0, rLX, 0, rLY, 0, rLZ);
    geometry.SetGradingCurve(NuTo::GeometryConcrete::B16, rNumClasses);
    geometry.SetParticleVolumeFraction(0.8);
    geometry.SetRelativeGrowthRate(0.1);
    geometry.SetInitialTimeBarrier(2);

    geometry.MaximizeParticleVolumeFraction(0.05);

    Eigen::MatrixXd particles = geometry.GetParticles(false);
    std::cout << "Created " << particles.rows() << " particles. " << std::endl;

    geometry.ExportGmshGeo3D(rGmshFile, 0.75);
}



int main(int argc, char* argv[])
{
    boost::filesystem::path outputPath = std::string(argv[0]) + "Out/";
    boost::filesystem::create_directory(outputPath);

    // When running with ctest, argv[1] contains the gmsh binary, something like '/usr/bin/gmsh'
    // However, 'gmsh' works fine in most cases.
    std::string gmshBinary = "gmsh";
    if (argc >= 2)
        gmshBinary = argv[1];

    std::cout << "Using gmsh binary: " << gmshBinary << std::endl;

    std::string gmshFile2D = outputPath.string() + "geometry2D";
    std::string gmshFile3D = outputPath.string() + "geometry3D";

    std::cout << "Gmsh File 2D:  " << gmshFile2D << ".geo" << std::endl;
    std::cout << "Gmsh File 3D:  " << gmshFile3D << ".geo" << std::endl;
    MaxDistanceMesh2D(gmshFile2D, 32, 16, 1);
    MaxDistanceMesh2D(gmshFile2D, 32, 16, 2);
    MaxDistanceMesh2D(gmshFile2D, 32, 16, 3);

    std::cout << "Meshing..." << std::endl;
    system((gmshBinary + " -3 -order 2 " + gmshFile2D + ".geo -o " + gmshFile2D + ".msh -v 2").c_str());

    MaxVolumeFraction3D(gmshFile3D, 32, 16, 16, 1);
    MaxVolumeFraction3D(gmshFile3D, 32, 16, 16, 2);
    MaxVolumeFraction3D(gmshFile3D, 32, 16, 16, 3);

}
