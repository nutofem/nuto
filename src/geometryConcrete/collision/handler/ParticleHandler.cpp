/*
 * ParticleHandler.cpp
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#include <fstream>
#include <iomanip>

#include "base/Exception.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"

#include "geometryConcrete/collision/handler/ParticleHandler.h"
#include "geometryConcrete/collision/collidables/CollidableParticleBase.h"
#include "geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "geometryConcrete/Specimen.h"

#include "geometryConcrete/GmshWriter.h"

#include "visualize/UnstructuredGrid.h"

NuTo::ParticleHandler::ParticleHandler(int rNumParticles, Eigen::MatrixXd rParticleBoundingBox, double rVelocityRange,
                                       double rGrowthRate, int rSeed)
    : mRNG(rSeed)
{
    mParticles.reserve(rNumParticles);

    double particleIndex = 100000;

    for (int i = 0; i < rNumParticles; ++i)
    {
        mParticles.push_back(new NuTo::CollidableParticleSphere(
                GetRandomVector(rParticleBoundingBox), GetRandomVector(-rVelocityRange / 2., rVelocityRange / 2.), 0.00,
                rGrowthRate, particleIndex++));
    }

    mVisualizationFileName = "spheres";
}

NuTo::ParticleHandler::ParticleHandler(Eigen::MatrixXd rSpheres, double rVelocityRange, double rRelativeGrowthRate,
                                       double rAbsoluteGrowthRate, int rSeed)
    : mRNG(rSeed)
{
    CreateParticlesFromMatrix(rSpheres, rVelocityRange, rRelativeGrowthRate, rAbsoluteGrowthRate);
    mVisualizationFileName = "spheres";
}

NuTo::ParticleHandler::ParticleHandler(const std::string& rFileName, double rVelocityRange, double rRelativeGrowthRate,
                                       double rAbsoluteGrowthRate, int rSeed)
    : mRNG(rSeed)
{
    Eigen::MatrixXd particles;
    NuTo::SerializeStreamIn sIn(rFileName, false);
    sIn >> particles;
    CreateParticlesFromMatrix(particles, rVelocityRange, rRelativeGrowthRate, rAbsoluteGrowthRate);
}

NuTo::ParticleHandler::~ParticleHandler()
{
    for (auto particle : mParticles)
        if (particle)
            delete particle;
}

Eigen::MatrixXd NuTo::ParticleHandler::GetParticles(bool rInitialRadius) const
{
    Eigen::MatrixXd particles(mParticles.size(), 4);
    int i = 0;
    for (auto particle : mParticles)
        particles.row(i++) = particle->ExportRow(rInitialRadius);

    return particles;
}

void NuTo::ParticleHandler::ExportParticlesToFile(const std::string& rExportFileName, bool rInitialRadius) const
{
    NuTo::SerializeStreamOut sOut(rExportFileName, false);
    sOut.SaveMatrix(GetParticles(rInitialRadius));
}

void NuTo::ParticleHandler::ExportParticlesToVTU3D(std::string rOutputDirectory, int rTimeStep, double rGlobalTime,
                                                   bool rFinal) const
{
    // modify rGlobalTime slightly to resolve every event
    rGlobalTime += rTimeStep * 1.e-10;

    NuTo::Visualize::UnstructuredGrid visuSpheres;
    visuSpheres.DefinePointData("Radius");
    visuSpheres.DefinePointData("Velocity");
    for (auto particle : mParticles)
        particle->VisualizationDynamic(visuSpheres, rFinal);

    std::stringstream fileName;
    fileName << rOutputDirectory << "/" << mVisualizationFileName << "_" << rTimeStep << ".vtu";
    visuSpheres.ExportVtuDataFile(fileName.str());

    // write an additional pvd file
    std::string resultFile = rOutputDirectory + std::string("/") + mVisualizationFileName + std::string(".pvd");
    std::fstream file;
    if (rTimeStep == 0)
        file.open(resultFile.c_str(), std::fstream::out);
    else
        file.open(resultFile.c_str(), std::fstream::out | std::fstream::in | std::ios_base::ate);

    if (!file.is_open())
        throw NuTo::Exception(std::string("[NuTo::CollisionHandler::VisualizeSpheres] Error opening file ") +
                              resultFile);

    std::stringstream endOfXML;
    endOfXML << "</Collection>" << std::endl;
    endOfXML << "</VTKFile>" << std::endl;
    if (rTimeStep == 0)
    {
        // header /////////////////////////////////////////////////////////////////
        file << "<?xml version=\"1.0\"?>" << std::endl;
        file << "<VTKFile type=\"Collection\">" << std::endl;
        file << "<Collection>" << std::endl;
    }
    else
    {
        // delete the last part of the xml file
        file.seekp(-endOfXML.str().length(), std::ios_base::end);
    }
    file << "<DataSet timestep=\"" << std::setprecision(10) << std::setw(10) << rGlobalTime << "\" file=\""
         << mVisualizationFileName << "_" << rTimeStep << ".vtu\"/>" << std::endl;

    file << endOfXML.str();
    file.close();
}

void NuTo::ParticleHandler::ExportParticlesToVTU2D(std::string rOutputFile, double rZCoord) const
{
    NuTo::Visualize::UnstructuredGrid visuSpheres;
    visuSpheres.DefinePointData("Radius");

    auto circles = GetParticles2D(rZCoord, 0);

    for (int i = 0; i < circles.rows(); i++)
    {
        Eigen::Vector3d coords;
        coords[0] = circles(i, 0);
        coords[1] = circles(i, 1);
        coords[2] = 0;

        unsigned int index = visuSpheres.AddPoint(coords);
        double radius = circles(i, 2);
        visuSpheres.SetPointData(index, "Radius", radius);
    }
    visuSpheres.ExportVtuDataFile(rOutputFile);
}

void NuTo::ParticleHandler::Sync(const double rTime)
{
    for (auto particle : mParticles)
        particle->MoveAndGrow(rTime);
}

double NuTo::ParticleHandler::GetKineticEnergy() const
{
    double eKin = 0.;
    for (auto particle : mParticles)
        eKin += particle->GetKineticEnergy();
    return eKin;
}

double NuTo::ParticleHandler::GetVolume() const
{
    double volume = 0.;
    for (auto particle : mParticles)
        volume += particle->GetVolume();
    return volume;
}

NuTo::CollidableParticleSphere* NuTo::ParticleHandler::GetParticle(const int rIndex) const
{
    return mParticles[rIndex];
}

int NuTo::ParticleHandler::GetNumParticles() const
{
    return mParticles.size();
}

double NuTo::ParticleHandler::GetAbsoluteMininimalDistance(Specimen& rSpecimen)
{

    struct StatBox
    {
        StatBox(Eigen::Matrix<double, 3, 2> rBox)
            : mBox(rBox){};
        std::vector<int> mSphereIndices;
        Eigen::Matrix<double, 3, 2> mBox;
    };

    Eigen::Vector3i divs = GetSubBoxDivisions(rSpecimen, 10);

    Eigen::Vector3d subBoxLength;
    for (int i = 0; i < 3; ++i)
        subBoxLength[i] = rSpecimen.GetLength(i) / divs[i];

    // create sub boxes
    std::vector<StatBox> boxes;
    for (int iX = 0; iX < divs[0]; ++iX)
        for (int iY = 0; iY < divs[1]; ++iY)
            for (int iZ = 0; iZ < divs[2]; ++iZ)
            {
                Eigen::Matrix<double, 3, 2> bounds;
                bounds << rSpecimen.GetBoundingBox()(0, 0) + iX * subBoxLength[0],
                        rSpecimen.GetBoundingBox()(0, 0) + (iX + 1) * subBoxLength[0],
                        rSpecimen.GetBoundingBox()(1, 0) + iY * subBoxLength[1],
                        rSpecimen.GetBoundingBox()(1, 0) + (iY + 1) * subBoxLength[1],
                        rSpecimen.GetBoundingBox()(2, 0) + iZ * subBoxLength[2],
                        rSpecimen.GetBoundingBox()(2, 0) + (iZ + 1) * subBoxLength[2];
                boxes.push_back(StatBox(bounds));
            }
    // add spheres

    for (unsigned int s = 0; s < mParticles.size(); ++s)
    {
        // get eight points to test, initialize with sphere centre
        std::vector<Eigen::Vector3d> corners(8, GetParticle(s)->GetPosition());
        // add radius to get different points
        int cornerIndex = 0;
        for (int iX = -1; iX <= 1; iX += 2)
            for (int iY = -1; iY <= 1; iY += 2)
                for (int iZ = -1; iZ <= 1; iZ += 2)
                {
                    double radius = GetParticle(s)->GetRadius0();
                    corners[cornerIndex](0) += iX * radius;
                    corners[cornerIndex](1) += iY * radius;
                    corners[cornerIndex](2) += iZ * radius;
                    cornerIndex++;
                }

        for (int c = 0; c < 8; ++c)
        {
            // get xyz cell index of corner points
            int indX = floor(static_cast<double>(corners[c](0) / subBoxLength[0]));
            int indY = floor(static_cast<double>(corners[c](1) / subBoxLength[1]));
            int indZ = floor(static_cast<double>(corners[c](2) / subBoxLength[2]));

            indX = std::min(indX, static_cast<int>(divs[0] - 1));
            indY = std::min(indY, static_cast<int>(divs[1] - 1));
            indZ = std::min(indZ, static_cast<int>(divs[2] - 1));

            indX = std::max(indX, 0);
            indY = std::max(indY, 0);
            indZ = std::max(indZ, 0);

            // add sphere to corresponding box
            int boxIndex = indX * divs[1] * divs[2] + indY * divs[2] + indZ;

            // only add new spheres
            if (boxes[boxIndex].mSphereIndices.size() != 0)
            {
                unsigned int lastSphereIndex = boxes[boxIndex].mSphereIndices.back();
                if (lastSphereIndex == s)
                    continue;
            }
            // sphere list is empty or lastIndex != boxIndex
            boxes[boxIndex].mSphereIndices.push_back(s);
        }
    }

    // only check boxes
    double minDistance = INFINITY;

    int nBoxes = divs[0] * divs[1] * divs[2];

    for (int b = 0; b < nBoxes; ++b)
    {
        StatBox& box = boxes[b];
        for (unsigned int i = 0; i < box.mSphereIndices.size(); ++i)
        {
            int indexI = box.mSphereIndices[i];
            for (unsigned int j = i + 1; j < box.mSphereIndices.size(); ++j)
            {
                int indexJ = box.mSphereIndices[j];

                Eigen::Vector3d dP = GetParticle(indexI)->GetPosition() - GetParticle(indexJ)->GetPosition();
                double dR = GetParticle(indexI)->GetRadius0() + GetParticle(indexJ)->GetRadius0();

                double ijDistance = sqrt(dP.dot(dP)) - dR;
                minDistance = std::min(minDistance, ijDistance);
            }
        }
    }
    return minDistance;
}

void NuTo::ParticleHandler::CreateParticlesFromMatrix(const Eigen::MatrixXd& rSpheres, double rVelocityRange,
                                                      double rRelativeGrowthRate, double rAbsoluteGrowthRate)
{
    int numRows = rSpheres.rows();
    mParticles.reserve(numRows);

    double particleIndex = 100000;

    for (int i = 0; i < numRows; ++i)
    {
        Eigen::Vector3d position({rSpheres(i, 0), rSpheres(i, 1), rSpheres(i, 2)});
        Eigen::Vector3d velocity = GetRandomVector(-rVelocityRange / 2., rVelocityRange / 2.);

        double radius = rSpheres(i, 3);
        mParticles.push_back(new NuTo::CollidableParticleSphere(
                position, velocity, radius, radius * rRelativeGrowthRate + rAbsoluteGrowthRate, particleIndex++));
    }
}

Eigen::VectorXd NuTo::ParticleHandler::GetRandomVector(double rStart, double rEnd)
{
    Eigen::Vector3d randomVector;
    std::uniform_real_distribution<double> distr(rStart, rEnd);

    for (int i = 0; i < 3; ++i)
        randomVector[i] = distr(mRNG);

    return randomVector;
}

void NuTo::ParticleHandler::ResetVelocities()
{
    for (auto particle : mParticles)
        particle->ResetVelocity();
}

Eigen::VectorXd NuTo::ParticleHandler::GetRandomVector(const Eigen::MatrixXd& rBounds)
{
    int vectorSize = rBounds.rows();
    Eigen::VectorXd randomVector(vectorSize);
    for (int i = 0; i < vectorSize; ++i)
    {
        std::uniform_real_distribution<double> distr(rBounds(i, 0), rBounds(i, 1));
        randomVector[i] = distr(mRNG);
    }
    return randomVector;
}

void NuTo::ParticleHandler::SetVisualizationFileName(const std::string& rVisualizationFileName)
{
    mVisualizationFileName = rVisualizationFileName;
}

Eigen::VectorXi NuTo::ParticleHandler::GetSubBoxDivisions(Specimen& rSpecimen, const int rParticlesPerBox)
{
    Eigen::VectorXi divs(3);

    double volBoundary = rSpecimen.GetVolume();
    double volSubBox = volBoundary / (GetNumParticles() / rParticlesPerBox);
    double lengthSubBox = std::pow(volSubBox, 1. / 3.);

    for (int i = 0; i < 3; ++i)
    {
        divs[i] = std::max(static_cast<int>(std::ceil(rSpecimen.GetLength(i) / lengthSubBox)), 1);
    }

    return divs;
}

Eigen::MatrixXd NuTo::ParticleHandler::GetParticles2D(double rZCoord, double rMinRadius) const
{
    Eigen::MatrixXd circles(1000, 3);
    int numCircles = 0;

    auto spheres = GetParticles(false);

    for (int countSphere = 0; countSphere < spheres.rows(); countSphere++)
    {
        double delta = spheres(countSphere, 2) - rZCoord;
        if (std::abs(delta) < spheres(countSphere, 3))
        {
            double radius =
                    sqrt(static_cast<double>(spheres(countSphere, 3) * spheres(countSphere, 3) - delta * delta));
            if (radius > rMinRadius)
            {
                // add circle
                if (numCircles == circles.rows())
                {
                    circles.conservativeResize(numCircles + 1000, 3);
                }
                circles(numCircles, 0) = spheres(countSphere, 0);
                circles(numCircles, 1) = spheres(countSphere, 1);
                circles(numCircles, 2) = radius;
                numCircles++;
            }
        }
    }
    circles.conservativeResize(numCircles, 3);

    return circles;
}

void NuTo::ParticleHandler::ExportParticlesToGmsh3D(std::string rOutputFile, Specimen& rSpecimen,
                                                    double rMeshSize) const
{
    GmshWriter::Options opts(rMeshSize, rMeshSize, 0);
    auto bounds = rSpecimen.GetBoundingBox();

    switch (rSpecimen.GetTypeOfSpecimen())
    {
    case Specimen::eSpecimenType::Cylinder:
    {
        double radius = bounds(0, 1);
        double height = bounds(2, 1);
        GmshWriter::Write(rOutputFile, GmshWriter::Cylinder{radius, height}, GetParticles(), opts);
        break;
    }
    case Specimen::eSpecimenType::Box:
    {
        GmshWriter::Box3D box(bounds.col(0), bounds.col(1));
        GmshWriter::Write(rOutputFile, box, GetParticles(), opts);
        break;
    }
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Don't know how to export this specimen type");
    }
}

void NuTo::ParticleHandler::ExportParticlesToGmsh2D(std::string rOutputFile, Specimen& rSpecimen, double rMeshSize,
                                                    double rZCoord, double rMinRadius) const
{
    GmshWriter::Options opts(rMeshSize, rMeshSize, 0);
    auto bounds = rSpecimen.GetBoundingBox();
    GmshWriter::Box2D box(bounds.block<2, 1>(0, 0), bounds.block<2, 1>(0, 1));
    GmshWriter::Write(rOutputFile, box, GetParticles2D(rZCoord, rMinRadius), opts);
}
