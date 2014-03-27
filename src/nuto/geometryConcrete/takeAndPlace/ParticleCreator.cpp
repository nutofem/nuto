/*
 * ParticleCreator.cpp
 *
 *  Created on: 27 Feb 2014
 *      Author: junger, improved by ttitsche
 */

#include "nuto/geometryConcrete/takeAndPlace/ParticleCreator.h"
#include "nuto/geometryConcrete/InputReader.h"
#include "nuto/base/Exception.h"
#include "nuto/geometryConcrete/WallTime.h"

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::ParticleCreator::CreateSpheresInSpecimen(
		const InputReader& rInput,
		const int rSeed)
{
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> spheresBoundary(0, 4);
	return CreateSpheresInSpecimen(
			rInput.GetBoxType(),
			rInput.GetBoundingBox(),
			rInput.GetVolumeFraction(),
			rInput.GetGradingCurve(),
			0.,
			rInput.GetAbsoluteDistance(),
			rSeed,
			spheresBoundary,
			rInput.GetShrinkage());
}


NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::ParticleCreator::CreateSpheresInSpecimen(
		const int rTypeOfSpecimen,
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox,
		const double rRelParticleVolume,
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rGradingCurve,
		const double rRelativeDistance,
		const double rAbsoluteDistance,
		const int rSeed,
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSpheresBoundary,
		const double rShrinkage)
{

	double tStart = WallTime::Get();

	// random number generator
	srand(rSeed);

	CheckBoundingBox(rBoundingBox);
	CheckGradingCurve(rGradingCurve);

	// calculate specimen length
	double volumeSpecimen = GetSpecimenVolume(rTypeOfSpecimen, rBoundingBox);

	FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> particles = PerformTakePhase(
			rGradingCurve,
			rSpheresBoundary,
			volumeSpecimen,
			rRelParticleVolume);

	double volume = 0.;
	double volumeShrinkage = 0.;

	for (int i = 0; i < particles.GetNumRows(); ++i)
	{
		volume += 4. / 3. * M_PI * particles(i, 3)*particles(i, 3)*particles(i, 3);
		particles(i, 3) = particles(i, 3) - rShrinkage;
		volumeShrinkage += 4. / 3. * M_PI * particles(i, 3)*particles(i, 3)*particles(i, 3);
	}


	std::cout << "[Take-Phase: ] Created " << particles.GetNumRows() << " particles. ";
	std::cout << "[Take-Phase: ] Phi = " << volume/volumeSpecimen << ", phi_shrinkage = " << volumeShrinkage/volumeSpecimen << std::endl;

	PerformPlacePhase(
			particles,
			rBoundingBox,
			rTypeOfSpecimen,
			rRelativeDistance,
			rAbsoluteDistance);


	std::cout << std::endl << "[Take-And-Place] Took " << WallTime::Get() - tStart << "s." << std::endl;

	return particles.GetBlock(rSpheresBoundary.GetNumRows(), 0, particles.GetNumRows() - rSpheresBoundary.GetNumRows(), 4);
}

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::ParticleCreator::PerformTakePhase(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rGradingCurve,
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSpheresBoundary,
		const double rVolumeSpecimen,
		const double rRelParticleVolume)
{
	// volume of particles per class
	int numGradingClasses = rGradingCurve.GetNumRows();

	std::vector<double> Vsoll(numGradingClasses);
	std::vector<double> Vist(numGradingClasses);
	std::vector<int> numParticlesPerClass(numGradingClasses);

	// calculating mass of the aggregates */
	double volumeSumParticles = rVolumeSpecimen * rRelParticleVolume;

	int numParticles = rSpheresBoundary.GetNumRows();
	FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> particles(0, 4);
	particles = rSpheresBoundary;


	for (int gc = 0; gc < numGradingClasses; ++gc)
	{
		double dMin = rGradingCurve(gc, 0);
		double dMax = rGradingCurve(gc, 1);
		double volumeFrac = rGradingCurve(gc, 2);

		numParticlesPerClass[gc] = 0;

		// calculate reference volume fraction of the mineral-size-class
		Vsoll[gc] = volumeSumParticles * volumeFrac;
		if (gc > 0)
		{
			Vsoll[gc] = Vsoll[gc] + (Vsoll[gc - 1] - Vist[gc - 1]);
		}
		Vist[gc] = 0.0;

		// generate particles until the reference volume fraction is reached
		bool finished = false;
		while (!finished)
		{
			// check size of table
			if (numParticles == particles.GetNumRows())
			{
				particles.ConservativeResizeRows(particles.GetNumRows() + 1000);
			}

			// calculate radius and volume of the particle
			double randomNumber = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // = geometry_rng.rand();

			double radius = 0.5 * dMin * dMax /
					pow((1.0 - randomNumber) * (dMax * dMax * dMax) + randomNumber * (dMin * dMin * dMin), 1.0 / 3.0);

			// volume
			double volumeParticle = 4.0 / 3.0 * M_PI * radius * radius * radius;

			Vist[gc] += volumeParticle;

			//create new particle
			if (Vist[gc] < Vsoll[gc])
			{
				//std::cout << "sphere " << numParticles+1 << " " << radius << " volume " << volumeParticle << std::endl;
				particles(numParticles, 3) = radius;
				numParticles++;
				numParticlesPerClass[gc]++;
			}
			else
			{
				finished = true;
				Vist[gc] -= volumeParticle;
			}
		}

		std::cout << "Volume for class " << gc + 1 << " : " << Vist[gc] / rVolumeSpecimen
				<< "(" << Vsoll[gc] / rVolumeSpecimen << ")" << std::endl;

	}

	particles.ConservativeResizeRows(numParticles);

	//sort
	std::sort(((double*) &particles.data()[3 * particles.GetNumRows()]),
			((double*) &particles.data()[3 * particles.GetNumRows() + numParticles]), std::greater<double>());

	return particles;

}

void NuTo::ParticleCreator::PerformPlacePhase(
		FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rParticles,
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox,
		const int rTypeOfSpecimen,
		const double rRelativeDistance,
		const double rAbsoluteDistance)
{

	int numParticles = rParticles.GetNumRows();
	FullVector<double, 3> boxLength = GetBoxLength(rBoundingBox);

	std::vector<double> sizeClasses = GetSizeClasses(rParticles);

	std::vector<double> numParticlesPerSize = GetNumParticlesPerSizeClass(rParticles, sizeClasses);

	int numParticlesAdded = 0;
	double lastPrintedFraction = 0.;
	for (unsigned int sizeClass = 0; sizeClass < sizeClasses.size(); ++sizeClass)
	{
//		std::cout << "Class D < " << sizeClasses[sizeClass] << " : " << numParticlesPerSize[sizeClass] << " particles." << std::endl;

		//create boxes for the previously inserted particles
		//width of each box = largest diameter
		FullVector<int, 3> nSubBox;
		FullVector<double, 3> lSubBox;
		for (int count = 0; count < 3; count++)
		{
			nSubBox[count] = std::floor(boxLength.GetValue(count) / sizeClasses[sizeClass]);
			lSubBox[count] = boxLength[count] / nSubBox[count];
		}

		int numberOfSubBoxes = nSubBox[0] * nSubBox[1] * nSubBox[2];
		std::vector<std::vector<int> > subBox(numberOfSubBoxes);

		// add all previously added particles to the new boxes
		for (int indexParticle = 0;
				indexParticle < numParticlesAdded;
				indexParticle++)
		{
			InsertParticleIntoBox(rParticles, indexParticle, subBox, nSubBox, lSubBox, rBoundingBox);
		}

		//loop all particles in the next new size class

		// now start inserting new particles into the box
		for (int countParticle = numParticlesAdded;
				countParticle < numParticlesAdded + numParticlesPerSize[sizeClass]; countParticle++)
		{
			bool inserted(false);
			int numTries(0);
			while (!inserted)
			{
				//create random coordinate
				FullVector<int, 3> cSubBox;
				for (int count = 0; count < 3; count++)
				{
					double radius = rParticles(countParticle, 3);
					rParticles(countParticle, count) = radius + (boxLength[count] - 2 * radius)
							* static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
					cSubBox[count] = (rParticles(countParticle, count) - rBoundingBox(count, 0)) / lSubBox[count];
				}

				//check for overlapping with the boundary
				if (CollidesWithBoundary(
						rBoundingBox,
						rTypeOfSpecimen,
						rParticles.GetRow(countParticle).transpose(),
						rRelativeDistance, rAbsoluteDistance))
					continue;

				//calculate the corresponding box
				int theBox = cSubBox[0] * nSubBox[1] * nSubBox[2] + cSubBox[1] * nSubBox[2] + cSubBox[2];

				//check for overlap with all the ellipses in that box
				bool noSeparation(true);
				for (unsigned int indexOther = 0;
						indexOther < subBox[theBox].size();
						indexOther++)
				{
					int other = subBox[theBox][indexOther];
					double deltaX = rParticles(countParticle, 0) - rParticles(other, 0);
					double deltaY = rParticles(countParticle, 1) - rParticles(other, 1);
					double deltaZ = rParticles(countParticle, 2) - rParticles(other, 2);
					double sumR = rParticles(countParticle, 3) * (1. + rRelativeDistance) +
							rAbsoluteDistance + rParticles(other, 3);
					if (deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ < sumR * sumR)
					{
						noSeparation = false;
						break;
					}
				}

				if (noSeparation)
				{
					//insert
					inserted = true;
					InsertParticleIntoBox(rParticles, countParticle, subBox, nSubBox, lSubBox, rBoundingBox);

					if ((double) (countParticle) / numParticles - lastPrintedFraction > 0.05)
					{
						if (lastPrintedFraction == 0.)
							std::cout << "[Take-Phase: ] Particles inserted: " << std::endl;

						std::cout << (double) (countParticle) / numParticles * 100.	<< "%..."<< std::endl;;
						lastPrintedFraction = (double) (countParticle) / numParticles;
					}

				}
				else
				{
					numTries++;
					if (numTries > 1e7)
						throw Exception("[NuTo::ParticleCreator::CreateSpheresInSpecimen] unable to insert sphere after a 1e7 tries.");
				}
			}
		}
		numParticlesAdded += numParticlesPerSize[sizeClass];
	}
}

void NuTo::ParticleCreator::CheckBoundingBox(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox)
		{
	if (rBoundingBox.GetNumRows() != 3 && rBoundingBox.GetNumColumns() != 2)
		throw Exception("[NuTo::ParticleCreator::CheckBoundingBox] bounding box has to have the dimension [3,2]");

	if (rBoundingBox.GetValue(0, 0) != 0. || rBoundingBox.GetValue(1, 0) != 0. || rBoundingBox.GetValue(2, 0) != 0.)
		throw Exception("[NuTo::ParticleCreator::CheckBoundingBox] bounding box has to start at (0.,0.,0.)");
}

void NuTo::ParticleCreator::CheckGradingCurve(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rGradingCurve)
		{
	int numGradingClasses = rGradingCurve.GetNumRows();

	if (numGradingClasses < 1)
		throw Exception(
				"[NuTo::ParticleCreator::CheckGradingCurve] at least one class in the grading curve should be defined.");

	for (int gc = 0; gc < numGradingClasses; ++gc)
	{
		double dMin = rGradingCurve(gc, 0);
		double dMax = rGradingCurve(gc, 1);
		if (dMin > dMax)
			throw Exception(
					"[NuTo::ParticleCreator::CheckGradingCurve] the minimum radius is larger than the maximum radius.");
		double volumeFrac = rGradingCurve(gc, 2);
		if (volumeFrac < 0 || volumeFrac > 1)
			throw Exception(
					"[NuTo::ParticleCreator::CheckGradingCurve] the mass fraction should be in the range [0,1].");
	}
}

NuTo::FullVector<double, 3> NuTo::ParticleCreator::GetBoxLength(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox)
		{
	FullVector<double, 3> lBox;
	for (int count = 0; count < 3; count++)
	{
		lBox[count] = rBoundingBox(count, 1) - rBoundingBox(count, 0);
		if (lBox[count] <= 0)
			throw Exception(
					"[NuTo::ParticleCreator::GetLBox] box dimensions should be not negative.");
	}
	return lBox;
}

const double NuTo::ParticleCreator::GetSpecimenVolume(
		int rTypeOfSpecimen,
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox)
		{
	double Vspecimen;

	FullVector<double, 3> lBox = GetBoxLength(rBoundingBox);


	switch (rTypeOfSpecimen)
	{
	case 0:
		//box
		Vspecimen = lBox[0] * lBox[1] * lBox[2];
		break;
	case 1:
		{
		Vspecimen = lBox[0] * lBox[1] * lBox[2];
		double D = rBoundingBox(0, 1) - rBoundingBox(0, 0);
		if (fabs(
				static_cast<double>(rBoundingBox(1, 1) - rBoundingBox(1, 0) - 1.5 * D)) > 1e-10)
			throw Exception("[NuTo::ParticleCreator::GetSpecimenVolume] "
					+ std::string("for the dog bone specimen, the y dimension should be 1.5 times the x dimension."));
		//subtract the circles
		double radius = 0.725 * D;
		double deltaAngle = 2. * 0.2 / 0.525;
		Vspecimen -= 2. * (deltaAngle / (2. * M_PI) * M_PI * radius * radius) * lBox[2];
	}
		break;
	case 2:
		{
		double D = rBoundingBox(0, 1) - rBoundingBox(0, 0);
		Vspecimen = M_PI * 0.25 * D * D * lBox[2];
		if (fabs(static_cast<double>(lBox[0] - lBox[1]) > 1e-10))
			throw Exception("[NuTo::ParticleCreator::GetSpecimenVolume] "
					+ std::string("for the cylindern, the x and y dimension should be identical (Diameter)."));
		if (D < 1e-10)
			throw Exception("[NuTo::ParticleCreator::GetSpecimenVolume] "
					+ std::string("for the cylindern, the x,y dimension should be positive (Diameter)."));
	}
		break;
	default:
		throw Exception("[NuTo::ParticleCreator::GetSpecimenVolume] specimen type not implemented.");
	}
	if (Vspecimen < 1.0e-14)
	{
		throw Exception("[NuTo::StructureBase::CreateSpheresInBox] negative volume of the box.");
	}
	return Vspecimen;
}


bool NuTo::ParticleCreator::CollidesWithBoundary(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox,
		const int rTypeOfSpecimen,
		const FullVector<double, 4> rParticle,
		const double rRelativeDistance,
		const double rAbsoluteDistance)
{
	bool collidesWithBoundary = false;
	switch (rTypeOfSpecimen)
	{
	case 0:
		break;
	case 1:
		{
		double D = rBoundingBox(0, 1) - rBoundingBox(0, 0);
		double deltaX = rParticle(0) - (rBoundingBox(0, 1) + 0.525 * D);
		double deltaY = rParticle(1) - (rBoundingBox(1, 0) + 0.75 * D);
		double sumR = rParticle(3) + 0.725 * D;
		if (deltaX * deltaX + deltaY * deltaY < sumR * sumR)
			collidesWithBoundary = true;

		deltaX = rParticle(0) - (rBoundingBox(0, 0) - 0.525 * D);
		deltaY = rParticle(1) - (rBoundingBox(1, 0) + 0.75 * D);
		if (deltaX * deltaX + deltaY * deltaY < sumR * sumR)
			collidesWithBoundary = true;
	}
		break;
	case 2:
		{
		double D = rBoundingBox(0, 1) - rBoundingBox(0, 0);
		double deltaX = rParticle(0) - (rBoundingBox(0, 0) + 0.5 * D);
		double deltaY = rParticle(1) - (rBoundingBox(1, 0) + 0.5 * D);
		double sumR = 0.5 * D - rParticle(3) * (1. + rRelativeDistance) - rAbsoluteDistance;
		if (sumR < 0)
			throw Exception("[NuTo::ParticleCreator::CreateSpheresInSpecimen] that should not have happend.");

		if (deltaX * deltaX + deltaY * deltaY > sumR * sumR)
			collidesWithBoundary = true;
	}
		break;
	default:
		throw Exception("[NuTo::ParticleCreator::CreateSpheresInSpecimen] specimen type not implemented.");
	}
	return collidesWithBoundary;
}

const std::vector<double> NuTo::ParticleCreator::GetSizeClasses(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rParticles)
{

	const double sizeScaleFactor = 1.5;
	double startDiameter = 2. * rParticles(rParticles.GetNumRows() - 1, 3);
	double endDiameter = 2. * rParticles(0, 3);
	double diameter = startDiameter;

	std::vector<double> sizesReverse;
	do
	{
		diameter *= sizeScaleFactor;
		sizesReverse.push_back(diameter);
	}while (diameter < endDiameter);

	// reverse vector
	std::vector<double> sizes(sizesReverse.size());
	for (unsigned int i = 0; i < sizesReverse.size(); ++i)
		sizes[sizes.size() - i - 1] = sizesReverse[i];

	return sizes;
}

const std::vector<double> NuTo::ParticleCreator::GetNumParticlesPerSizeClass(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rParticles,
		const std::vector<double>& rSizes)
{

	unsigned int numParticles = rParticles.GetNumRows();

	std::vector<double> numParticlesPerSize(rSizes.size());

	unsigned int sizeCount = 0;
	unsigned int sizeClass = 0;
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		double diameter = 2 * rParticles(i, 3);
		if (diameter >= rSizes[sizeClass + 1])
			sizeCount++;
		else
		{
			numParticlesPerSize[sizeClass] = sizeCount;
			sizeClass++;
			sizeCount = 1;
			if (sizeClass == rSizes.size() - 1)
			{
				numParticlesPerSize[sizeClass] = numParticles - i;
				break;
			}
		}
	}
	return numParticlesPerSize;
}


//! @brief cut spheres at a given z-coordinate to create circles (in 2D)
//! @parameters rSpheres matrix with the spheres (x,y,z,r)
//! @parameters rZCoord z coordinate (where to cut)
//! @parameters rMinRadius minimal radius of the circle
//! @return ... matrix with the circles (x,y,r)
NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::ParticleCreator::CutSpheresZ(
		NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSpheres,
		double rZCoord, double rMinRadius)
{
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> circles(1000, 3);
	int numCircles(0);
	for (int countSphere = 0; countSphere < rSpheres.GetNumRows();
			countSphere++)
	{
		double delta = rSpheres(countSphere, 2) - rZCoord;
		if (fabs(delta) < rSpheres(countSphere, 3))
		{
			double radius = sqrt(static_cast<double>(rSpheres(countSphere, 3) * rSpheres(countSphere, 3) - delta * delta));
			if (radius > rMinRadius)
			{
				//add circle
				if (numCircles == circles.GetNumRows())
				{
					circles.ConservativeResizeRows(numCircles + 1000);
				}
				circles(numCircles, 0) = rSpheres(countSphere, 0);
				circles(numCircles, 1) = rSpheres(countSphere, 1);
				circles(numCircles, 2) = radius;
				numCircles++;
			}
		}
	}
	circles.ConservativeResizeRows(numCircles);

	return circles;
}

//! @brief ... inserts a particle into subboxes to increase efficiency when performing overlap checks
void NuTo::ParticleCreator::InsertParticleIntoBox(
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rParticles,
		const int rTheParticle,
		std::vector<std::vector<int> >& rSubBox,
		const FullVector<int, 3>& rNSubBox,
		const FullVector<double, 3>& rLSubBox,
		const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox)
		{

	FullVector<int, 3> cSubBoxMin;
	FullVector<int, 3> cSubBoxMax;
	for (int coordinate = 0; coordinate < 3; coordinate++)
	{
		double radius = rParticles(rTheParticle, 3);
		double posParticle = rParticles(rTheParticle, coordinate);
		double posBoundary = rBoundingBox(coordinate, 0);


		// bigger particles may extend one box
		// --> calculate the index range in each direction
		int indexMin = (posParticle - radius - posBoundary) / rLSubBox[coordinate];
		int indexMax = (posParticle + radius - posBoundary) / rLSubBox[coordinate];

		// the diameter of newly added particles is at least rLSubBox[coordinate]
		// As it may lay on a box border, the index range of existing boxes
		// has to be extended by 1

		indexMin--;
		indexMax++;

		// prevent segmentation faults
		cSubBoxMin[coordinate] = std::max(indexMin, 0);
		cSubBoxMax[coordinate] = std::min(indexMax, rNSubBox.GetValue(coordinate) - 1);

	}

	//insert the center box + all the surroundings
	for (int iX = cSubBoxMin[0]; iX <= cSubBoxMax[0]; iX++)
		for (int iY = cSubBoxMin[1]; iY <= cSubBoxMax[1]; iY++)
			for (int iZ = cSubBoxMin[2]; iZ <= cSubBoxMax[2]; iZ++)
			{
				int theBox = iX * rNSubBox[1] * rNSubBox[2] + iY * rNSubBox[2] + iZ;
				rSubBox[theBox].push_back(rTheParticle);
			}
}

