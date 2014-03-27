/*
 * ParticleHandler.cpp
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/collision/handler/ParticleHandler.h"
#include "nuto/math/FullVector.h"


NuTo::ParticleHandler::ParticleHandler(ParticleHandler& rOther)
		: mParticleIndex(rOther.mParticleIndex)
{
	mParticles = rOther.mParticles;
}

NuTo::ParticleHandler::ParticleHandler(
		const int rNumParticles,
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBoundingBox,
		const double rVelocityRange,
		const double rGrowthRate)
{
	mParticleIndex = 100000;

	mParticles.resize(rNumParticles);

	for (int i = 0; i < rNumParticles; ++i)
	{
		mParticles[i] = new NuTo::CollidableParticleSphere(
				GetRandomVector(rBoundingBox),
				GetRandomVector(-rVelocityRange / 2, rVelocityRange / 2),
				0.,
				rGrowthRate,
				mParticleIndex++);
	}

}

NuTo::ParticleHandler::ParticleHandler(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rSpheres,
		const double rVelocityRange,
		const double rGrowthRate)
{
	mParticleIndex = 100000;

	int numRows = rSpheres.GetNumRows();

	mParticles.resize(numRows);

	for (int i = 0; i < numRows; ++i)
	{
		NuTo::FullVector<double, 3> position( { rSpheres(i, 0), rSpheres(i, 1), rSpheres(i, 2) });
		NuTo::FullVector<double, 3> velocity = GetRandomVector(-rVelocityRange / 2, rVelocityRange / 2);
		double radius = rSpheres(i, 3);
		mParticles[i] = new NuTo::CollidableParticleSphere(
				position,
				velocity,
				radius,
				rGrowthRate,
				mParticleIndex++);
	}
}

NuTo::ParticleHandler::~ParticleHandler()
{

}

void NuTo::ParticleHandler::Delete()
{
	for (auto particle : mParticles)
		delete particle;

}

void NuTo::ParticleHandler::AddParticle(const CollidableParticleSphere& rSphere)
{
	mParticles.push_back(new CollidableParticleSphere(rSphere));
}

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::ParticleHandler::GetParticles() const
{
	FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> particles(mParticles.size(), 4);
	int i = 0;
	for (auto particle : mParticles)
		particles.SetRow(i++, particle->ExportRow());

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
	fileName << rName << "/spheres_" << rTimeStep << ".vtu";
	visuSpheres.ExportVtuDataFile(fileName.str());

	//write an additional pvd file
	std::string resultFile = rName + std::string("/spheres.pvd");
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
			<< rGlobalTime << "\" file=\"spheres_" << rTimeStep << ".vtu\"/>"
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


double NuTo::ParticleHandler::GetAbsoluteMininimalDistance(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBoundingBox)
{
	double rMax = GetRMax();

	FullVector<double, 3> length;
	FullVector<double, 3> subLength;
	FullVector<int, 3> divs;
	for (int i = 0; i < 3; ++i)
	{
		length[i] = rBoundingBox(i, 1) - rBoundingBox(i, 0);
		divs[i] = static_cast<int>(length[i] / (rMax));
		subLength[i] = length[i] / divs[i];
	}

	// create sub boxes
	std::vector<StatBox> boxes;
	for (int iX = 0; iX < divs[0]; ++iX)
		for (int iY = 0; iY < divs[1]; ++iY)
			for (int iZ = 0; iZ < divs[2]; ++iZ)
			{
				FullMatrix<double, 3, 2> bounds;
				bounds << rBoundingBox(0, 0) + iX * subLength[0],
						rBoundingBox(0, 0) + (iX + 1) * subLength[0],
						rBoundingBox(1, 0) + iY * subLength[1],
						rBoundingBox(1, 0) + (iY + 1) * subLength[1],
						rBoundingBox(2, 0) + iZ * subLength[2],
						rBoundingBox(2, 0) + (iZ + 1) * subLength[2];
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
			int indX = floor(static_cast<double>(corners[c](0) / subLength[0]));
			int indY = floor(static_cast<double>(corners[c](1) / subLength[1]));
			int indZ = floor(static_cast<double>(corners[c](2) / subLength[2]));

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
				double ijDistance = dP.Dot(dP) - dR * dR;
				minDistance = std::min(minDistance, ijDistance);
			}
		}
	}
	return sqrt(minDistance);

}

double NuTo::ParticleHandler::GetRMax()
{
	double rMax = 0;
	for (auto particle : mParticles)
		rMax = std::max(rMax, particle->GetRadius0());
	return rMax;
}

NuTo::FullVector<double, 3> NuTo::ParticleHandler::GetRandomVector(const double rStart, const double rEnd)
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

NuTo::FullVector<double, 3> NuTo::ParticleHandler::GetRandomVector(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBounds)
{
	NuTo::FullVector<double, 3> randomVector;
	for (int i = 0; i < 3; ++i)
	{
		double rnd = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
		rnd *= rBounds.GetValue(i,1) - rBounds.GetValue(i,0);   // correct range
		rnd += rBounds.GetValue(i,0);			// correct starting point
		randomVector.SetValue(i, rnd);
	}
	return randomVector;
}

