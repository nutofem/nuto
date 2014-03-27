/*
 * ParticleCreator.h
 *
 *  Created on: 27 Feb 2014
 *      Author: ttitsche
 */

#ifndef PARTICLECREATOR_H_
#define PARTICLECREATOR_H_

#include "nuto/math/FullMatrix.h"

namespace NuTo
{
class InputReader;
class ParticleCreator
{
public:

	//! @brief ... creates randomly distributed, non-overlapping particles
	//! @param rTypeOfSpecimen ... 0 box, 1 dogbone, 2 cylinder
	//! @param rBoundingBox ... box for the spheres (3*2 matrix)
	//! @param rRelParticleVolume ... percentage of particle volume inside the box
	//! @param rGradingCurve ... matrix with each line min_diameter, max_diameter, volume percentage of that sieve size
	//! @param rRelativeDistance ... scaling factor to increase the diameter when inserting the sphere to ensure a minimum distance
	//! @param rAbsoluteDistance ... scaling value to increase the diameter when inserting the sphere to ensure a minimum distance
	//! @param rSeed ... seed for the random number generator
	//! @param rSpheresBoundary ... particles simulated on the boundary e.g. created with CreateSpheresOnBoxBoundary (they do not contribute to the grading curve)
	//! @param rShrinkage ... absolute value (mm) of particle shrinkage, allows temporary higher particle densities
	//! @return ... matrix with spheres (coordinates x y z and radius)
	static FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> CreateSpheresInSpecimen(
			const int rTypeOfSpecimen,
			const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBoundingBox,
			const double rRelParticleVolume,
			const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rGradingCurve,
			const double rRelativeDistance,
			const double rAbsoluteDistance,
			const int rSeed,
			const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rSpheresBoundary,
			const double rShrinkage);

	//! @brief ... creates randomly distributed, non-overlapping particles
	//! @param rInput ... NuTo::InputReader object
	//! @param rSeed ... seed for the random number generator
	//! @return ... matrix with spheres (coordinates x y z and radius)
	static FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> CreateSpheresInSpecimen(
			const InputReader& rInput, const int rSeed);



	//! @brief ... cut spheres at a given z-coordinate to create circles (in 2D)
	//! @param rSpheres matrix with the spheres (x,y,z,r)
	//! @param rZCoord z coordinate (where to cut)
	//! @param rMinRadius minimal radius of the circle
	//! @return ... matrix with the circles (x,y,r)
	static FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> CutSpheresZ(
			FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSpheres,
			double rZCoord, double rMinRadius);

	//! @brief ... calculates the volume of the specimen
	//! @param rTypeOfSpecimen ... 0 box, 1 dogbone, 2 cylinder
	//! @param rBoundingBox ... box for the spheres (3*2 matrix)
	//! @return volume of the specimen
	static const double GetSpecimenVolume(
			int rTypeOfSpecimen,
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox);

	//! @brief ... performs the "take"-of the "take-and-place" algorithm
	//! @return ... matrix with sphere radii according to the grading curve
	static FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> PerformTakePhase(
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rGradingCurve,
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSpheresBoundary,
			const double rVolumeSpecimen,
			const double rRelParticleVolume);

	//! @brief ... performs the "place"-of the "take-and-place" algorithm
	//! @return ... matrix with randomly distributed, non-overlapping particles
	static void PerformPlacePhase(
			FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rParticles,
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox,
			const int rTypeOfSpecimen,
			const double rRelativeDistance,
			const double rAbsoluteDistance);

private:



	//! @brief ... inserts a particle into subboxes to increase efficiency when performing overlap checks
	static void InsertParticleIntoBox(
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rParticles,
			const int rTheParticle,
			std::vector<std::vector<int> >& rSubBox,
			const FullVector<int, 3>& rNSubBox,
			const FullVector<double, 3>& rLSubBox,
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox);

	static void CheckBoundingBox(
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox);

	static void CheckGradingCurve(
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rGradingCurve);

	static FullVector<double, 3> GetBoxLength(
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox);

	//! @brief ... collision check with specimen boundary
	static bool CollidesWithBoundary(
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBoundingBox,
			const int rTypeOfSpecimen,
			const FullVector<double, 4> rParticle,
			const double rRelativeDistance,
			const double rAbsoluteDistance);

	//! @brief ... recalculate size classes != grading curve classes for better performance
	static const std::vector<double> GetSizeClasses(
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rParticles);

	static const std::vector<double> GetNumParticlesPerSizeClass(
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rParticles,
			const std::vector<double>& rSizes);
};

} /* namespace NuTo */
#endif /* PARTICLECREATOR_H_ */
