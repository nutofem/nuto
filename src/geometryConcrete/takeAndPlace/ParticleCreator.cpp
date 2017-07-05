/*
 * ParticleCreator.cpp
 *
 *  Created on: 27 Feb 2014
 *      Author: junger, improved by ttitsche
 */
#include "base/Exception.h"
#include "base/Timer.h"

#include "geometryConcrete/takeAndPlace/ParticleCreator.h"
#include "geometryConcrete/InputReader.h"


NuTo::ParticleCreator::ParticleCreator(NuTo::Specimen rSpecimen, const double rShrinkage, const long rNumMaxTries)
    : mSpecimen(rSpecimen)
    , mShrinkage(rShrinkage)
    , mNumMaxTries(rNumMaxTries)
{
    mVolume = mSpecimen.GetVolume();
}

Eigen::MatrixXd NuTo::ParticleCreator::CreateSpheresInSpecimen(const double rRelParticleVolume,
                                                               const Eigen::MatrixXd& rGradingCurve,
                                                               const double rRelativeDistance,
                                                               const double rAbsoluteDistance, const int rSeed,
                                                               const Eigen::MatrixXd& rSpheresBoundary) const
{
    Timer timer(__FUNCTION__, true);

    // random number generator
    srand(rSeed);

    CheckGradingCurve(rGradingCurve);

    Eigen::MatrixXd particles = PerformTakePhase(rGradingCurve, rSpheresBoundary, rRelParticleVolume);

    PerformPlacePhase(particles, rRelativeDistance, rAbsoluteDistance);

    return particles.block(rSpheresBoundary.rows(), 0, particles.rows() - rSpheresBoundary.rows(), 4);
}

Eigen::MatrixXd NuTo::ParticleCreator::PerformTakePhase(const Eigen::MatrixXd& rGradingCurve,
                                                        const Eigen::MatrixXd& rSpheresBoundary,
                                                        const double rRelParticleVolume) const
{
    // volume of particles per class
    int numGradingClasses = rGradingCurve.rows();

    std::vector<double> Vsoll(numGradingClasses);
    std::vector<double> Vist(numGradingClasses);
    std::vector<int> numParticlesPerClass(numGradingClasses);

    // calculating mass of the aggregates */
    double volumeSumParticles = mVolume * rRelParticleVolume;

    int numParticles = rSpheresBoundary.rows();
    Eigen::MatrixXd particles(0, 4);
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
            if (numParticles == particles.rows())
            {
                particles.conservativeResize(particles.rows() + 1000, 4);
            }

            // calculate radius and volume of the particle
            double randomNumber = static_cast<double>(rand()) / static_cast<double>(RAND_MAX); // = geometry_rng.rand();

            double radius =
                    0.5 * dMin * dMax /
                    pow((1.0 - randomNumber) * (dMax * dMax * dMax) + randomNumber * (dMin * dMin * dMin), 1.0 / 3.0);

            // volume
            double volumeParticle = 4.0 / 3.0 * M_PI * radius * radius * radius;

            Vist[gc] += volumeParticle;

            // create new particle
            if (Vist[gc] < Vsoll[gc])
            {
                // std::cout << "sphere " << numParticles+1 << " " << radius << " volume " << volumeParticle <<
                // std::endl;
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

        std::cout << "Volume for class " << gc + 1 << " : " << Vist[gc] / mVolume << "(" << Vsoll[gc] / mVolume << ")"
                  << numParticles << " particles." << std::endl;
    }

    particles.conservativeResize(numParticles, 4);

    // sort
    std::sort(((double*)&particles.data()[3 * particles.rows()]),
              ((double*)&particles.data()[3 * particles.rows() + numParticles]), std::greater<double>());

    double volume = 0.;
    double volumeShrinkage = 0.;

    for (int i = 0; i < particles.rows(); ++i)
    {
        volume += GetVolume(particles(i, 3));
        particles(i, 3) = particles(i, 3) * (1 - mShrinkage);
        volumeShrinkage += GetVolume(particles(i, 3));
    }

    std::cout << "[Take-Phase: ] Created " << particles.rows() << " particles. ";
    std::cout << "[Take-Phase: ] Phi = " << volume / mVolume << ", phi_shrinkage = " << volumeShrinkage / mVolume
              << std::endl;

    return particles;
}

void NuTo::ParticleCreator::PerformPlacePhase(Eigen::MatrixXd& rParticles, const double rRelativeDistance,
                                              const double rAbsoluteDistance) const
{

    auto numParticles = rParticles.rows();

    std::vector<double> sizeClasses = GetSizeClasses(rParticles);

    std::vector<double> numParticlesPerSize = GetNumParticlesPerSizeClass(rParticles, sizeClasses);


    int numParticlesAdded = 0;
    double lastPrintedFraction = 0.;
    for (unsigned int sizeClass = 0; sizeClass < sizeClasses.size(); ++sizeClass)
    {
        std::cout << "Class D < " << sizeClasses[sizeClass] << " : " << numParticlesPerSize[sizeClass] << " particles."
                  << std::endl;

        // create boxes for the previously inserted particles
        // width of each box = largest diameter
        Eigen::Vector3i nSubBox;
        Eigen::Vector3d lSubBox;
        for (int count = 0; count < 3; count++)
        {
            // rescale to original sizes
            double subBoxLength = sizeClasses[sizeClass] / (1. - mShrinkage);
            nSubBox[count] = std::floor(mSpecimen.GetLength()[count] / subBoxLength);
            lSubBox[count] = mSpecimen.GetLength()[count] / nSubBox[count];
        }

        int numberOfSubBoxes = nSubBox[0] * nSubBox[1] * nSubBox[2];
        std::vector<std::vector<int>> subBox(numberOfSubBoxes);

        // add all previously added particles to the new boxes
        for (int indexParticle = 0; indexParticle < numParticlesAdded; indexParticle++)
        {
            InsertParticleIntoBox(rParticles, indexParticle, subBox, nSubBox, lSubBox);
        }

        // loop all particles in the next new size class

        // now start inserting new particles into the box
        for (int countParticle = numParticlesAdded; countParticle < numParticlesAdded + numParticlesPerSize[sizeClass];
             countParticle++)
        {
            bool inserted = false;
            int numTries = 0;
            while (!inserted)
            {
                // create random coordinate
                Eigen::Vector3i cSubBox;
                for (int count = 0; count < 3; count++)
                {
                    double radius = rParticles(countParticle, 3);
                    double randomNumber = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

                    rParticles(countParticle, count) = radius + mSpecimen.GetBoundingBox()(count, 0) +
                                                       (mSpecimen.GetLength()[count] - 2. * radius) * randomNumber;
                    cSubBox[count] =
                            (rParticles(countParticle, count) - mSpecimen.GetBoundingBox()(count, 0)) / lSubBox[count];
                }

                //				check for overlapping with the boundary
                if (not mSpecimen.IsBox())
                    if (CollidesWithBoundary(rParticles.row(countParticle).transpose(), rRelativeDistance,
                                             rAbsoluteDistance))
                        continue;

                // calculate the corresponding box
                int theBox = cSubBox[0] * nSubBox[1] * nSubBox[2] + cSubBox[1] * nSubBox[2] + cSubBox[2];

                // check for overlap with all the ellipses in that box
                bool noSeparation(true);
                for (int other : subBox[theBox])
                {
                    double deltaX = rParticles(countParticle, 0) - rParticles(other, 0);
                    double deltaY = rParticles(countParticle, 1) - rParticles(other, 1);
                    double deltaZ = rParticles(countParticle, 2) - rParticles(other, 2);
                    double sumR = rParticles(countParticle, 3) * (1. + rRelativeDistance) + rAbsoluteDistance +
                                  rParticles(other, 3);
                    if (deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ < sumR * sumR)
                    {
                        noSeparation = false;
                        break;
                    }
                }

                if (noSeparation)
                {
                    // insert
                    inserted = true;
                    InsertParticleIntoBox(rParticles, countParticle, subBox, nSubBox, lSubBox);

                    if ((double)(countParticle) / numParticles - lastPrintedFraction > 0.05)
                    {
                        if (lastPrintedFraction == 0.)
                            std::cout << "[Take-Phase: ] Particles inserted: " << std::endl;

                        std::cout << (double)(countParticle) / numParticles * 100. << "%..." << std::endl;
                        lastPrintedFraction = (double)(countParticle) / numParticles;
                    }
                }
                else
                {
                    numTries++;
                    if (numTries > mNumMaxTries)
                    {
                        std::stringstream exceptionStream;
                        exceptionStream
                                << "[NuTo::ParticleCreator::CreateSpheresInSpecimen] unable to insert sphere after ";
                        exceptionStream << mNumMaxTries << " tries.";
                        throw Exception(exceptionStream.str());
                    }
                }
            }
        }
        numParticlesAdded += numParticlesPerSize[sizeClass];
    }
}

void NuTo::ParticleCreator::CheckGradingCurve(const Eigen::MatrixXd& rGradingCurve) const
{
    auto numGradingClasses = rGradingCurve.rows();

    if (numGradingClasses < 1)
        throw Exception("[NuTo::ParticleCreator::CheckGradingCurve] at least one class in the grading curve should be "
                        "defined.");

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


bool NuTo::ParticleCreator::CollidesWithBoundary(const Eigen::Vector4d& rParticle, const double rRelativeDistance,
                                                 const double rAbsoluteDistance) const
{

    auto& bBox = mSpecimen.GetBoundingBox();

    bool collidesWithBoundary = false;
    switch (mSpecimen.GetTypeOfSpecimen())
    {
    case Specimen::Box:
        break;
    case Specimen::Dogbone:
    {
        double D = bBox(0, 1) - bBox(0, 0);
        double deltaX = rParticle(0) - (bBox(0, 1) + 0.525 * D);
        double deltaY = rParticle(1) - (bBox(1, 0) + 0.75 * D);
        double sumR = rParticle(3) + 0.725 * D;
        if (deltaX * deltaX + deltaY * deltaY < sumR * sumR)
            collidesWithBoundary = true;

        deltaX = rParticle(0) - (bBox(0, 0) - 0.525 * D);
        deltaY = rParticle(1) - (bBox(1, 0) + 0.75 * D);
        if (deltaX * deltaX + deltaY * deltaY < sumR * sumR)
            collidesWithBoundary = true;
    }
    break;
    case Specimen::Cylinder:
    {
        double D = bBox(0, 1) - bBox(0, 0);
        double deltaX = rParticle(0) - (bBox(0, 0) + 0.5 * D);
        double deltaY = rParticle(1) - (bBox(1, 0) + 0.5 * D);
        double sumR = 0.5 * D - rParticle(3) * (1. + rRelativeDistance) - rAbsoluteDistance;
        if (sumR < 0)
            throw Exception(__PRETTY_FUNCTION__, "that should not have happend.");

        if (deltaX * deltaX + deltaY * deltaY > sumR * sumR)
            collidesWithBoundary = true;
    }
    break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "specimen type not implemented.");
    }
    return collidesWithBoundary;
}

const std::vector<double> NuTo::ParticleCreator::GetSizeClasses(const Eigen::MatrixXd& rParticles) const
{

    const double sizeScaleFactor = 2;
    double startDiameter = 2. * rParticles(rParticles.rows() - 1, 3);
    double endDiameter = 2. * rParticles(0, 3);
    double diameter = startDiameter;

    std::vector<double> sizesReverse;
    do
    {
        diameter *= sizeScaleFactor;
        sizesReverse.push_back(diameter);
    } while (diameter < endDiameter);

    // reverse vector
    std::vector<double> sizes(sizesReverse.size());
    for (unsigned int i = 0; i < sizesReverse.size(); ++i)
        sizes[sizes.size() - i - 1] = sizesReverse[i];

    std::cout << sizes.size() << std::endl;

    // reset biggest size class to biggest particle diameter
    sizes[0] = 2. * rParticles(0, 3);

    return sizes;
}

const std::vector<double> NuTo::ParticleCreator::GetNumParticlesPerSizeClass(const Eigen::MatrixXd& rParticles,
                                                                             const std::vector<double>& rSizes) const
{

    unsigned int numParticles = rParticles.rows();

    std::vector<double> numParticlesPerSize(rSizes.size());

    if (rSizes.size() == 1)
    {
        numParticlesPerSize[0] = numParticles;
        return numParticlesPerSize;
    }

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

//! @brief ... inserts a particle into subboxes to increase efficiency when performing overlap checks
void NuTo::ParticleCreator::InsertParticleIntoBox(const Eigen::MatrixXd& rParticles, const int rTheParticle,
                                                  std::vector<std::vector<int>>& rSubBox,
                                                  const Eigen::Vector3i& rNSubBox,
                                                  const Eigen::Vector3d& rLSubBox) const
{

    Eigen::Vector3i cSubBoxMin;
    Eigen::Vector3i cSubBoxMax;

    int coordMax = 3;

    for (int coordinate = 0; coordinate < coordMax; coordinate++)
    {
        double radius = rParticles(rTheParticle, 3);
        double posParticle = rParticles(rTheParticle, coordinate);
        double posBoundary = mSpecimen.GetBoundingBox()(coordinate, 0);

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
        cSubBoxMax[coordinate] = std::min(indexMax, rNSubBox[coordinate] - 1);
    }

    // insert the center box + all the surroundings
    for (int iX = cSubBoxMin[0]; iX <= cSubBoxMax[0]; iX++)
        for (int iY = cSubBoxMin[1]; iY <= cSubBoxMax[1]; iY++)
            for (int iZ = cSubBoxMin[2]; iZ <= cSubBoxMax[2]; iZ++)
            {
                int theBox = iX * rNSubBox[1] * rNSubBox[2] + iY * rNSubBox[2] + iZ;
                rSubBox[theBox].push_back(rTheParticle);
            }
}
double NuTo::ParticleCreator::GetVolume(double radius) const
{
    return 4. / 3. * M_PI * radius * radius * radius;
}
