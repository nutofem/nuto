/*
 * ParticleCreator.h
 *
 *  Created on: 27 Feb 2014
 *      Author: ttitsche
 */

#pragma once


#include <vector>
#include "geometryConcrete/Specimen.h"


namespace NuTo
{
class InputReader;

class ParticleCreator
{
public:

	//! @param rSpecimen ... Specimen object
	ParticleCreator(NuTo::Specimen rSpecimen, const double rShrinkage, const long rNumMaxTries = 10000000);

	//! @brief ... creates randomly distributed, non-overlapping particles
	//! @param rRelParticleVolume ... percentage of particle volume inside the box
	//! @param rGradingCurve ... matrix with each line min_diameter, max_diameter, volume percentage of that sieve size
	//! @param rRelativeDistance ... scaling factor to increase the diameter when inserting the sphere to ensure a minimum distance
	//! @param rAbsoluteDistance ... scaling value to increase the diameter when inserting the sphere to ensure a minimum distance
	//! @param rSeed ... seed for the random number generator
	//! @param rSpheresBoundary ... particles simulated on the boundary e.g. created with CreateSpheresOnBoxBoundary (they do not contribute to the grading curve)
	//! @param rShrinkage ... absolute value (mm) of particle shrinkage, allows temporary higher particle densities
	//! @return ... matrix with spheres (coordinates x y z and radius)
	Eigen::MatrixXd CreateSpheresInSpecimen(
			const double rRelParticleVolume,
			const Eigen::MatrixXd& rGradingCurve,
			const double rRelativeDistance,
			const double rAbsoluteDistance,
			const int rSeed,
			const Eigen::MatrixXd& rSpheresBoundary) const;

	//! @brief ... performs the "take"-of the "take-and-place" algorithm
	//! @return ... matrix with sphere radii according to the grading curve
	Eigen::MatrixXd PerformTakePhase(
			const Eigen::MatrixXd& rGradingCurve,
			const Eigen::MatrixXd& rSpheresBoundary,
			const double rRelParticleVolume) const;

	//! @brief ... performs the "place"-of the "take-and-place" algorithm
	//! @return ... matrix with randomly distributed, non-overlapping particles
	void PerformPlacePhase(
			Eigen::MatrixXd& rParticles,
			const double rRelativeDistance,
			const double rAbsoluteDistance) const;

private:

	//! @brief ... inserts a particle into subboxes to increase efficiency when performing overlap checks
	void InsertParticleIntoBox(
			const Eigen::MatrixXd& rParticles,
			const int rTheParticle,
			std::vector<std::vector<int> >& rSubBox,
			const Eigen::Vector3i& rNSubBox,
			const Eigen::Vector3d& rLSubBox) const;

	void CheckGradingCurve(
			const Eigen::MatrixXd& rGradingCurve) const;

	//! @brief ... collision check with specimen boundary
	bool CollidesWithBoundary(
			const Eigen::Vector4d& rParticle,
			const double rRelativeDistance,
			const double rAbsoluteDistance) const;

	//! @brief ... recalculate size classes != grading curve classes for better performance
	const std::vector<double> GetSizeClasses(
			const Eigen::MatrixXd& rParticles) const;

	const std::vector<double> GetNumParticlesPerSizeClass(
			const Eigen::MatrixXd& rParticles,
			const std::vector<double>& rSizes) const;
	double GetVolume(double radius) const;

	NuTo::Specimen mSpecimen;
	const double mShrinkage;
	const long mNumMaxTries;

	double mVolume;
};

} /* namespace NuTo */
