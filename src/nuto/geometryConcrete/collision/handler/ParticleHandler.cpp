/*
 * ParticleHandler.cpp
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/collision/handler/ParticleHandler.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleBase.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "nuto/geometryConcrete/Specimen.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#endif

NuTo::ParticleHandler::ParticleHandler(
		const int rNumParticles,
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rParticleBoundingBox,
		const double rVelocityRange,
		const double rGrowthRate)
		:
				mParticleIndex(100000)
{
	mParticles.reserve(rNumParticles);

	for (int i = 0; i < rNumParticles; ++i)
	{
		mParticles.push_back(new NuTo::CollidableParticleSphere(
				GetRandomVector(rParticleBoundingBox),
				GetRandomVector(-rVelocityRange / 2., rVelocityRange / 2.),
				0.00,
				rGrowthRate,
				mParticleIndex++));
	}

	mVisualizationFileName = "spheres";
}

NuTo::ParticleHandler::ParticleHandler(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rSpheres,
		const double rVelocityRange,
		const double rRelativeGrowthRate,
		const double rAbsoluteGrowthRate)
		:
				mParticleIndex(100000)
{

	int numRows = rSpheres.GetNumRows();

	mParticles.reserve(numRows);

	bool is2D = Is2DSimulation(rSpheres);
	std::cout << is2D << std::endl;

	for (int i = 0; i < numRows; ++i)
	{
		NuTo::FullVector<double, 3> position( { rSpheres(i, 0), rSpheres(i, 1), rSpheres(i, 2) });
		NuTo::FullVector<double, 3> velocity = GetRandomVector(-rVelocityRange / 2., rVelocityRange / 2.);
		if (is2D)
			velocity[2] = 0.;
		double radius = rSpheres(i, 3);
		mParticles.push_back(new NuTo::CollidableParticleSphere(
				position,
				velocity,
				radius,
				radius * rRelativeGrowthRate + rAbsoluteGrowthRate,
				mParticleIndex++));
	}
	mVisualizationFileName = "spheres";
}

NuTo::ParticleHandler::~ParticleHandler()
{
	for (auto particle : mParticles)
		if (particle)
			delete particle;
}

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::ParticleHandler::GetParticles(bool rInitialRadius) const
{
	FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> particles(mParticles.size(), 4);
	int i = 0;
	for (auto particle : mParticles)
		particles.SetRow(i++, particle->ExportRow(rInitialRadius));

	return particles;
}

void NuTo::ParticleHandler::VisualizeSpheres(std::string rName, int rTimeStep, double rGlobalTime, bool rFinal)
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
	fileName << rName << "/" << mVisualizationFileName << "_" << rTimeStep << ".vtu";
	visuSpheres.ExportVtuDataFile(fileName.str());

	//write an additional pvd file
	std::string resultFile = rName + std::string("/") + mVisualizationFileName + std::string(".pvd");
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

NuTo::CollidableParticleSphere* NuTo::ParticleHandler::GetParticle(const int rIndex)
{
	return mParticles[rIndex];
}

const int NuTo::ParticleHandler::GetNumParticles() const
{
	return mParticles.size();
}

double NuTo::ParticleHandler::GetAbsoluteMininimalDistance(Specimen& rSpecimen)
{

	struct StatBox
	{
		StatBox(FullMatrix<double, 3, 2> rBox)
				: mBox(rBox)
		{
		}
		;
		std::vector<int> mSphereIndices;
		FullMatrix<double, 3, 2> mBox;
	};

	FullVector<int, 3> divs = GetSubBoxDivisions(rSpecimen, 10);

	FullVector<double, 3> subBoxLength;
	for (int i = 0; i < 3; ++i)
		subBoxLength[i] = rSpecimen.GetLength(i) / divs[i];

	// create sub boxes
	std::vector<StatBox> boxes;
	for (int iX = 0; iX < divs[0]; ++iX)
		for (int iY = 0; iY < divs[1]; ++iY)
			for (int iZ = 0; iZ < divs[2]; ++iZ)
			{
				FullMatrix<double, 3, 2> bounds;
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
		std::vector<FullVector<double, 3> > corners(8, GetParticle(s)->GetPosition());
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

				FullVector<double, 3> dP = GetParticle(indexI)->GetPosition() - GetParticle(indexJ)->GetPosition();
				double dR = GetParticle(indexI)->GetRadius0() + GetParticle(indexJ)->GetRadius0();

				double ijDistance = sqrt(dP.Dot(dP)) - dR;
				minDistance = std::min(minDistance, ijDistance);
			}
		}
	}
	return minDistance;

}

NuTo::FullVector<double, Eigen::Dynamic> NuTo::ParticleHandler::GetRandomVector(const double rStart, const double rEnd)
{
	NuTo::FullVector<double, 3> randomVector;
	for (int i = 0; i < 3; ++i)
	{
		double rnd = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
		rnd *= rEnd - rStart;   // correct range
		rnd += rStart;			// correct starting point
		randomVector.SetValue(i, rnd);
	}
	return randomVector;
}

void NuTo::ParticleHandler::ResetVelocities()
{
	for (auto particle : mParticles)
		particle->ResetVelocity();
}

NuTo::FullVector<double, Eigen::Dynamic> NuTo::ParticleHandler::GetRandomVector(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBounds)
{
	int vectorSize = rBounds.GetNumRows();
	NuTo::FullVector<double, Eigen::Dynamic> randomVector(vectorSize);
	for (int i = 0; i < vectorSize; ++i)
	{
		double rnd = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
		rnd *= rBounds.GetValue(i, 1) - rBounds.GetValue(i, 0);   // correct range
		rnd += rBounds.GetValue(i, 0);			// correct starting point
		randomVector.SetValue(i, rnd);
	}
	return randomVector;
}

void NuTo::ParticleHandler::SetVisualizationFileName(const std::string& visualizationFileName)
{
	mVisualizationFileName = visualizationFileName;
}

NuTo::FullVector<int, Eigen::Dynamic> NuTo::ParticleHandler::GetSubBoxDivisions(Specimen& rSpecimen, const int rParticlesPerBox)
{

	FullVector<int, Eigen::Dynamic> divs(3);

	if (rSpecimen.Is2D())
	{
		double areaBoundary = rSpecimen.GetVolume();
		double areaSubBox = areaBoundary / (GetNumParticles() / rParticlesPerBox);
		double lengthSubBox = std::sqrt(areaSubBox);
		divs[0] = rSpecimen.GetLength(0) / lengthSubBox;
		divs[1] = rSpecimen.GetLength(1) / lengthSubBox;
		divs[2] = 1;
	}
	else
	{
		double volBoundary = rSpecimen.GetVolume();
		double volSubBox = volBoundary / (GetNumParticles() / rParticlesPerBox);
		double lengthSubBox = std::pow(volSubBox, 1. / 3.);

		for (int i = 0; i < 3; ++i)
			divs.SetValue(i, std::ceil(rSpecimen.GetLength(i) / lengthSubBox));
	}

	return divs;

}

bool NuTo::ParticleHandler::Is2DSimulation(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rSpheres)
{
	double coordZ = rSpheres(0, 2);
	for (int i = 0; i < rSpheres.GetNumRows(); ++i)
		if (coordZ != rSpheres(i, 2))
			return false;
	return true;
}
