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



#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeUnstructuredGrid.h"
#endif

NuTo::ParticleHandler::ParticleHandler(
		int rNumParticles,
		Eigen::MatrixXd rParticleBoundingBox,
		double rVelocityRange,
		double rGrowthRate,
        int rSeed)
    : mRNG(rSeed)
{
	mParticles.reserve(rNumParticles);

    double particleIndex = 100000;

	for (int i = 0; i < rNumParticles; ++i)
	{
		mParticles.push_back(new NuTo::CollidableParticleSphere(
				GetRandomVector(rParticleBoundingBox),
				GetRandomVector(-rVelocityRange / 2., rVelocityRange / 2.),
				0.00,
				rGrowthRate,
                particleIndex++));
	}

	mVisualizationFileName = "spheres";
}

NuTo::ParticleHandler::ParticleHandler(
		Eigen::MatrixXd rSpheres,
		double rVelocityRange,
		double rRelativeGrowthRate,
        double rAbsoluteGrowthRate,
        int rSeed)
    : mRNG(rSeed)
{
    CreateParticlesFromMatrix(rSpheres, rVelocityRange, rRelativeGrowthRate, rAbsoluteGrowthRate);
    mVisualizationFileName = "spheres";
}

NuTo::ParticleHandler::ParticleHandler(
    const std::string &rFileName,
    double rVelocityRange,
    double rRelativeGrowthRate,
    double rAbsoluteGrowthRate,
    int rSeed)
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

void NuTo::ParticleHandler::ExportParticlesToFile(const std::string &rExportFileName, bool rInitialRadius) const
{
    NuTo::SerializeStreamOut sOut(rExportFileName, false);
    sOut.SaveMatrix(GetParticles(rInitialRadius));
}

void NuTo::ParticleHandler::ExportParticlesToVTU3D(std::string rOutputDirectory, int rTimeStep, double rGlobalTime, bool rFinal) const
{
#ifdef ENABLE_VISUALIZE
	// modify rGlobalTime slightly to resolve every event
	rGlobalTime += rTimeStep * 1.e-10;

	NuTo::VisualizeUnstructuredGrid visuSpheres;
	visuSpheres.DefinePointDataScalar("Radius");
	visuSpheres.DefinePointDataVector("Velocity");
	for (auto particle : mParticles)
		particle->VisualizationDynamic(visuSpheres, rFinal);

	std::stringstream fileName;
    fileName << rOutputDirectory << "/" << mVisualizationFileName << "_" << rTimeStep << ".vtu";
	visuSpheres.ExportVtuDataFile(fileName.str());

	//write an additional pvd file
    std::string resultFile = rOutputDirectory + std::string("/") + mVisualizationFileName + std::string(".pvd");
	std::fstream file;
	if (rTimeStep == 0)
		file.open(resultFile.c_str(), std::fstream::out);
	else
		file.open(resultFile.c_str(), std::fstream::out | std::fstream::in | std::ios_base::ate);

	if (!file.is_open())
		throw NuTo::Exception(std::string("[NuTo::CollisionHandler::VisualizeSpheres] Error opening file ") + resultFile);

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
		//delete the last part of the xml file
		file.seekp(-endOfXML.str().length(), std::ios_base::end);
	}
	file << "<DataSet timestep=\"" << std::setprecision(10) << std::setw(10)
            << rGlobalTime << "\" file=\"" << mVisualizationFileName << "_" << rTimeStep << ".vtu\"/>"
			<< std::endl;

	file << endOfXML.str();
	file.close();
#endif
}

void NuTo::ParticleHandler::ExportParticlesToVTU2D(std::string rOutputFile, double rZCoord) const
{
#ifdef ENABLE_VISUALIZE

    NuTo::VisualizeUnstructuredGrid visuSpheres;
    visuSpheres.DefinePointDataScalar("Radius");

    auto circles = GetParticles2D(rZCoord,0);

    for (int i = 0; i < circles.rows(); i ++)
    {
        double coords[3];
        coords[0] = circles(i,0);
        coords[1] = circles(i,1);
        coords[2] = 0;

        unsigned int index = visuSpheres.AddPoint(coords);
        double radius = circles(i,2);
        visuSpheres.SetPointDataScalar(index, "Radius", radius);
    }
    visuSpheres.ExportVtuDataFile(rOutputFile);

#endif
}

void NuTo::ParticleHandler::Sync(const double rTime)
{
	for (auto particle : mParticles)
		particle->MoveAndGrow(rTime);

}

const double NuTo::ParticleHandler::GetKineticEnergy() const
{
	double eKin = 0.;
	for (auto particle : mParticles)
		eKin += particle->GetKineticEnergy();
	return eKin;

}

const double NuTo::ParticleHandler::GetVolume() const
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
				: mBox(rBox)
		{
		}
		;
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
				bounds <<
						rSpecimen.GetBoundingBox()(0, 0) + iX * subBoxLength[0],
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

void NuTo::ParticleHandler::CreateParticlesFromMatrix(const Eigen::MatrixXd& rSpheres,
                                                      double rVelocityRange,
                                                      double rRelativeGrowthRate,
                                                      double rAbsoluteGrowthRate)
{
    int numRows = rSpheres.rows();
    mParticles.reserve(numRows);

    double particleIndex = 100000;

    for (int i = 0; i < numRows; ++i)
    {
        Eigen::Vector3d position( { rSpheres(i, 0), rSpheres(i, 1), rSpheres(i, 2) });
        Eigen::Vector3d velocity = GetRandomVector(-rVelocityRange / 2., rVelocityRange / 2.);

        double radius = rSpheres(i, 3);
        mParticles.push_back(new NuTo::CollidableParticleSphere(
                position,
                velocity,
                radius,
                radius * rRelativeGrowthRate + rAbsoluteGrowthRate,
                particleIndex++));
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

Eigen::VectorXd NuTo::ParticleHandler::GetRandomVector(
		const Eigen::MatrixXd& rBounds)
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

Eigen::MatrixXd NuTo::ParticleHandler::GetParticles2D(
        double rZCoord, double rMinRadius) const
        {
    Eigen::MatrixXd circles(1000, 3);
    int numCircles = 0;

    auto spheres = GetParticles(false);

    for (int countSphere = 0; countSphere < spheres.rows();
            countSphere++)
    {
        double delta = spheres(countSphere, 2) - rZCoord;
        if (std::abs(delta) < spheres(countSphere, 3))
        {
            double radius = sqrt(static_cast<double>(spheres(countSphere, 3) * spheres(countSphere, 3) - delta * delta));
            if (radius > rMinRadius)
            {
                //add circle
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
    circles.conservativeResize(numCircles,3);

    return circles;
}

void NuTo::ParticleHandler::ExportParticlesToGmsh3D(std::string rOutputFile,
        Specimen& rSpecimen, double rMeshSize) const
{
    std::ofstream mFile;

    // Define header
    // some default options, adjust this afterwards in the .geo file

    mFile.open(rOutputFile);
    mFile << "Mesh.Algorithm = 6;\n";
    //mFile << "Mesh.HighOrderPoissonRatio = 0.2;\n";
    mFile << "Mesh.HighOrderOptimize = 1;\n"; //number of smoothing steps
    //mFile << "Mesh.HighOrderNumLayers = 6;\n";
    //mFile << "Mesh.HighOrderSmoothingThreshold = 0.2;\n"; //default 0.5
    //mFile << "Mesh.MultiplePassesMeshes = 1;\n";//default 0 do a simple mesh and use for background meshing
    mFile << "Mesh.Optimize = 2;\n";
    //mFile << "Mesh.CharacteristicLengthFromCurvatur = 0;\n";
    //mFile << "Mesh.OptimizeNetgen = 1;\n";
    //mFile << "Mesh.SecondOrderExperimental = 1;\n"; //experimental code for second order
    //mFile << "Mesh.SecondOrderLinear = 1;\n"; //should second order elements just be created by linear interpolation?
    //mFile << "Mesh.ThirdOrderLinear = " << true << ";\n"; //should second order elements just be created by linear interpolation?
    mFile << "Mesh.Smoothing = 2;\n"; //number of smoothing steps for the final mesh


    // define the sphere function

    mFile << "Function MySphere \n";
    mFile << "  meshCircleR = meshCircle; \n";
    mFile << "  p1 = newp; Point(p1) = {xC,  yC,  zC,  meshCircleR} ; \n";
    mFile << "  p2 = newp; Point(p2) = {xC+R,yC,  zC,  meshCircleR} ; \n";
    mFile << "  p3 = newp; Point(p3) = {xC,  yC+R,zC,  meshCircleR} ; \n";
    mFile << "  p4 = newp; Point(p4) = {xC,  yC,  zC+R,meshCircleR} ; \n";
    mFile << "  p5 = newp; Point(p5) = {xC-R,yC,  zC,  meshCircleR} ; \n";
    mFile << "  p6 = newp; Point(p6) = {xC,  yC-R,zC,  meshCircleR} ; \n";
    mFile << "  p7 = newp; Point(p7) = {xC,  yC,  zC-R,meshCircleR} ; \n";
    mFile << " \n";
    mFile << "  c1  = newreg; Circle(c1)  = {p2,p1,p7}; c2  = newreg; Circle(c2)  = {p7,p1,p5}; \n";
    mFile << "  c3  = newreg; Circle(c3)  = {p5,p1,p4}; c4  = newreg; Circle(c4)  = {p4,p1,p2}; \n";
    mFile << "  c5  = newreg; Circle(c5)  = {p2,p1,p3}; c6  = newreg; Circle(c6)  = {p3,p1,p5}; \n";
    mFile << "  c7  = newreg; Circle(c7)  = {p5,p1,p6}; c8  = newreg; Circle(c8)  = {p6,p1,p2}; \n";
    mFile << "  c9  = newreg; Circle(c9)  = {p7,p1,p3}; c10 = newreg; Circle(c10) = {p3,p1,p4}; \n";
    mFile << "  c11 = newreg; Circle(c11) = {p4,p1,p6}; c12 = newreg; Circle(c12) = {p6,p1,p7}; \n";
    mFile << " \n";
    mFile << "  l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1}; \n";
    mFile << "  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2}; \n";
    mFile << "  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3}; \n";
    mFile << "  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4}; \n";
    mFile << "  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5}; \n";
    mFile << "  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6}; \n";
    mFile << "  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7}; \n";
    mFile << "  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8}; \n";
    mFile << "   \n";
    mFile << "  theLoops[t] = newreg; \n";
    mFile << " \n";
    mFile << "  Surface Loop(theLoops[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1}; \n";
    mFile << " \n";
    mFile << "  thehole = newreg; \n";
    mFile << "  Volume(thehole) = theLoops[t]; \n";
    mFile << "  theAggregates[t] = thehole; \n";
    mFile << " \n";
    mFile << "Return \n";
    mFile << " \n";

    // define the specimen

    auto bounds = rSpecimen.GetBoundingBox();

    mFile << "xS = " << bounds(0, 0) << "; xE = " << bounds(0, 1) << ";\n";
    mFile << "yS = " << bounds(1, 0) << "; yE = " << bounds(1, 1) << ";\n";
    mFile << "zS = " << bounds(2, 0) << "; zE = " << bounds(2, 1) << ";\n";
    mFile << "meshSpecimen = " << rMeshSize << "; \n";

    mFile << "// defines a box-shaped specimen \n";
    mFile << "// by start coordinates <xyz>S \n";
    mFile << "// and end coordinates  <xyz>E \n";
    mFile << "\n";
    mFile << "// top points: z = zE \n";
    mFile << "p0 = newp; Point(p0) = {xS, yS, zE, meshSpecimen}; \n";
    mFile << "p1 = newp; Point(p1) = {xS, yE, zE, meshSpecimen}; \n";
    mFile << "p2 = newp; Point(p2) = {xE, yE, zE, meshSpecimen}; \n";
    mFile << "p3 = newp; Point(p3) = {xE, yS, zE, meshSpecimen}; \n";
    mFile << "// bottom points z = zS \n";
    mFile << "p4 = newp; Point(p4) = {xS, yS, zS, meshSpecimen}; \n";
    mFile << "p5 = newp; Point(p5) = {xS, yE, zS, meshSpecimen}; \n";
    mFile << "p6 = newp; Point(p6) = {xE, yE, zS, meshSpecimen}; \n";
    mFile << "p7 = newp; Point(p7) = {xE, yS, zS, meshSpecimen}; \n";
    mFile << "\n";
    mFile << "// top lines z = zS \n";
    mFile << "lT0 = newreg; Line(lT0) = {p0, p1}; \n";
    mFile << "lT1 = newreg; Line(lT1) = {p1, p2}; \n";
    mFile << "lT2 = newreg; Line(lT2) = {p2, p3}; \n";
    mFile << "lT3 = newreg; Line(lT3) = {p3, p0}; \n";
    mFile << "// bottom lines z = zE \n";
    mFile << "lB0 = newreg; Line(lB0) = {p4, p5}; \n";
    mFile << "lB1 = newreg; Line(lB1) = {p5, p6}; \n";
    mFile << "lB2 = newreg; Line(lB2) = {p6, p7}; \n";
    mFile << "lB3 = newreg; Line(lB3) = {p7, p4}; \n";
    mFile << "// connection zS --> zE \n";
    mFile << "lC0 = newreg; Line(lC0) = {p0, p4}; \n";
    mFile << "lC1 = newreg; Line(lC1) = {p1, p5}; \n";
    mFile << "lC2 = newreg; Line(lC2) = {p2, p6}; \n";
    mFile << "lC3 = newreg; Line(lC3) = {p3, p7}; \n";
    mFile << "\n";
    mFile << "// lineloops and surfaces \n";
    mFile << "// the index says which coordinate is constant \n";
    mFile << "lxS = newreg; Line Loop(lxS) = {-lT3, lC3, lB3,-lC0}; Plane Surface(newreg) = {lxS}; \n";
    mFile << "lxE = newreg; Line Loop(lxE) = {-lT1, lC1, lB1,-lC2}; Plane Surface(newreg) = {lxE}; \n";
    mFile << "\n";
    mFile << "lyS = newreg; Line Loop(lyS) = {-lT0, lC0, lB0,-lC1}; Plane Surface(newreg) = {lyS}; \n";
    mFile << "lyE = newreg; Line Loop(lyE) = {-lT2, lC2, lB2,-lC3}; Plane Surface(newreg) = {lyE}; \n";
    mFile << "\n";
    mFile << "lzS = newreg; Line Loop(lzS) = { lT0, lT1, lT2, lT3}; Plane Surface(newreg) = {lzS}; \n";
    mFile << "lzE = newreg; Line Loop(lzE) = {-lB3,-lB2,-lB1,-lB0}; Plane Surface(newreg) = {lzE}; \n";
    mFile << "\n";
    mFile << "theLoops[0] = newreg; \n";
    mFile << "Surface Loop(theLoops[0]) = {lxS+1, lyS+1, lxE+1, lyE+1, lzS+1, lzE+1}; \n";
    mFile << "theBox = newreg; \n";
    mFile << "Volume(theBox) = theLoops[0]; \n";
    mFile << "\n";
    mFile << "\n";

    auto spheres = GetParticles();
    int objectCounter = 1; // The first (index 0) object is the box.

    mFile << "meshCircle = " << rMeshSize << "; \n";
    for (int i = 0; i < spheres.rows(); i++)
    {
        mFile << "t = " << objectCounter << ";\n";
        objectCounter++;
        mFile << "xC = " << spheres(i, 0) << "; yC = " << spheres(i, 1)
                << "; zC = " << spheres(i, 2) << ";\n";
        mFile << "R = " << spheres(i, 3) << "; \n";
        mFile << "Call MySphere; \n";
        mFile << " \n";
        mFile << " \n";
    }


    mFile << "volNr = newreg; \n";
    mFile << "Volume(volNr) = {theLoops[]};\n";
    mFile << "Physical Volume(newreg) = volNr;\n";
    mFile << "Physical Volume(newreg) = {theAggregates[]};\n";
    mFile.close();
}

void NuTo::ParticleHandler::ExportParticlesToGmsh2D(std::string rOutputFile,
        Specimen& rSpecimen, double rMeshSize, double rZCoord,
        double rMinRadius) const
{
    std::ofstream file;

    // Define header
    // some default options, adjust this afterwards in the .geo file

    file.open(rOutputFile);
    file << "Mesh.Algorithm = 1;\n"; //Frontal hex
    //mFile << "Mesh.HighOrderPoissonRatio = 0.2;\n";
    //mFile << "Mesh.HighOrderOptimize = 1;\n"; //number of smoothing steps
    //mFile << "Mesh.HighOrderNumLayers = 6;\n";
    //mFile << "Mesh.HighOrderSmoothingThreshold = 0.2;\n"; //default 0.5
    //mFile << "Mesh.MultiplePassesMeshes = 1;\n";//default 0 do a simple mesh and use for background meshing
    file << "Mesh.Optimize = 1;\n";
    //mFile << "Mesh.CharacteristicLengthFromCurvatur = 0;\n";
    //mFile << "Mesh.OptimizeNetgen = 1;\n";
    //mFile << "Mesh.SecondOrderExperimental = 1;\n"; //experimental code for second order
    //mFile << "Mesh.SecondOrderLinear = 1;\n"; //should second order elements just be created by linear interpolation?
    //mFile << "Mesh.ThirdOrderLinear = " << true << ";\n"; //should second order elements just be created by linear interpolation?
    file << "Mesh.Smoothing = 2;\n"; //number of smoothing steps for the final mesh


    // define the circle function

    file << "Function MySphere \n";
    file << "  meshCircleR = meshCircle; \n";
    file << "  p1 = newp; Point(p1) = {xC,  yC,  0,  meshCircleR} ; \n";
    file << "  p2 = newp; Point(p2) = {xC+R,yC,  0,  meshCircleR} ; \n";
    file << "  p3 = newp; Point(p3) = {xC-R,yC,  0,  meshCircleR} ; \n";
    file << " \n";
    file << "  c1  = newreg; Circle(c1)  = {p2,p1,p3}; c2  = newreg; Circle(c2)  = {p3,p1,p2}; \n";
    file << " \n";
    file << "  l1 = newreg; Line Loop(l1) = {c1,c2};\n";
    file << "  s1 = newreg; Plane Surface(s1) = {l1}; \n";
    file << "   \n";
    file << "  theLoops[t] = l1; \n";
    file << "  theAggregates[t] = s1; \n";
    file << " \n";
    file << "Return \n";
    file << " \n";

    // define the specimen

    auto bounds = rSpecimen.GetBoundingBox();

    file << "xS = " << bounds(0, 0) << "; xE = " << bounds(0, 1) << ";\n";
    file << "yS = " << bounds(1, 0) << "; yE = " << bounds(1, 1) << ";\n";
    file << "meshSpecimen = " << rMeshSize << ";\n";
    file << "// defines a box-shaped specimen \n";
    file << "// by start coordinates <xyz>S \n";
    file << "// and end coordinates  <xyz>E \n";
    file << "\n";
    file << "// points: \n";
    file << "p0 = newp; Point(p0) = {xS, yS, 0, meshSpecimen}; \n";
    file << "p1 = newp; Point(p1) = {xE, yS, 0, meshSpecimen}; \n";
    file << "p2 = newp; Point(p2) = {xE, yE, 0, meshSpecimen}; \n";
    file << "p3 = newp; Point(p3) = {xS, yE, 0, meshSpecimen}; \n";
    file << "\n";
    file << "// lines \n";
    file << "l0 = newreg; Line(l0) = {p0, p1}; \n";
    file << "l1 = newreg; Line(l1) = {p1, p2}; \n";
    file << "l2 = newreg; Line(l2) = {p2, p3}; \n";
    file << "l3 = newreg; Line(l3) = {p3, p0}; \n";
    file << "\n";
    file << "// lineloops and surfaces \n";
    file << "// the index says which coordinate is constant \n";
    file << "box = newreg; Line Loop(box) = { l0, l1, l2, l3}; \n";
    file << "\n";
    file << "\n";

    // call the circle function

    auto circles = GetParticles2D(rZCoord, rMinRadius);

    if (circles.rows() == 0)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Found no aggregates for visualization. Change the z-Slice or increase rMin.");

    int objectCounter = 0;

    file << "meshCircle = " << rMeshSize << "; \n";

    for (int i = 0; i < circles.rows(); i++)
    {
        file << "t = " << objectCounter << ";\n";
        objectCounter++;
        file << "xC = " << circles(i, 0) << "; yC = " << circles(i, 1) << ";\n";
        file << "R = " << circles(i, 2) << "; \n";
        file << "Call MySphere; \n";
        file << " \n";
        file << " \n";

    }

    // finish

    file << "volNr = newreg; \n";
    file << "Plane Surface(volNr) = {box, theLoops[]};\n";
    file << "Physical Surface(newreg) = volNr;\n";
    file << "Physical Surface(newreg) = {theAggregates[]}; \n";
    file.close();
}
