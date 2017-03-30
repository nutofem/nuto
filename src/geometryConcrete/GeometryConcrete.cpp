/*
 * GeometryConcrete.cpp
 *
 *  Created on: 4 Sep 2015
 *      Author: ttitsche
 */

#include <iostream>

#include "base/Exception.h"

#include "geometryConcrete/GeometryConcrete.h"
#include "geometryConcrete/Specimen.h"
#include "geometryConcrete/takeAndPlace/ParticleCreator.h"
#include "geometryConcrete/collision/handler/ParticleHandler.h"
#include "geometryConcrete/collision/handler/SubBoxHandler.h"
#include "geometryConcrete/collision/handler/CollisionHandler.h"



NuTo::GeometryConcrete::GeometryConcrete()
{
}

NuTo::GeometryConcrete::~GeometryConcrete()
{
    delete mParticleHandler;
    delete mSpecimen;
}

void NuTo::GeometryConcrete::MaximizeParticleDistance(double rParticleDistance)
{

    CheckParameters();

    if (mAbsoluteGrowthRate <= 0)
        throw NuTo::Exception("[NuTo::GeometryConcrete::ExportGmsh2D] Call SetAbsoluteGrowthRate(rate) with rate > 0 first!");

    NuTo::ParticleCreator creator(*mSpecimen, 0);
    Eigen::MatrixXd spheresBoundary(0, 4);
    Eigen::MatrixXd spheres =
            creator.CreateSpheresInSpecimen(mParticleVolumeFraction, mGradingCurve, 0., 0., mSeed, spheresBoundary);

    SetParticles(spheres);

    NuTo::SubBoxHandler subBoxes(*mParticleHandler, *mSpecimen, 1);
    NuTo::CollisionHandler handler(*mParticleHandler, subBoxes, "");


    double timeEnd = .5*rParticleDistance/mAbsoluteGrowthRate;

    try
    {
        handler.Simulate(mNumEventsMax, timeEnd, mSecondsWallTimeMax, mSecondsPrint, mInitialTimeBarrier);
    } catch (NuTo::Exception& e)
    {
        e.AddMessage("The simulation stopped with an exception. \n");
        if (mContinueOnException)
        {
            std::cout << e.ErrorMessage() << "\n but I'll continue.";
        }
        else
        {
            throw;
        }
    }

    std::cout << "min Dist.= " << mParticleHandler->GetAbsoluteMininimalDistance(*mSpecimen)<< std::endl;

    // the spheres are stored with the diameters _before_ the EDMD simulation
    NuTo::ParticleHandler* old = mParticleHandler;
    mParticleHandler = new ParticleHandler(old->GetParticles(true), mRandomVelocityRange, mRelativeGrowthRate, mAbsoluteGrowthRate);

    delete old;
}

void NuTo::GeometryConcrete::MaximizeParticleVolumeFraction(double rShrinkage)
{

    CheckParameters();

    if (mRelativeGrowthRate <= 0)
        throw NuTo::Exception("[NuTo::GeometryConcrete::MaximizeParticleVolumeFraction] Call SetRelativeGrowthRate(rate) with rate > 0 first!");

    NuTo::ParticleCreator creator(*mSpecimen, rShrinkage);
    Eigen::MatrixXd spheresBoundary(0, 4);
    Eigen::MatrixXd spheres =
            creator.CreateSpheresInSpecimen(mParticleVolumeFraction, mGradingCurve, 0., 0., mSeed, spheresBoundary);

    SetParticles(spheres);

    NuTo::SubBoxHandler subBoxes(*mParticleHandler, *mSpecimen, 10);
    NuTo::CollisionHandler handler(*mParticleHandler, subBoxes, "");

    double timeEnd = (1. / (1. - rShrinkage) - 1.) / mRelativeGrowthRate;

    try
    {
        handler.Simulate(mNumEventsMax, timeEnd, mSecondsWallTimeMax, mSecondsPrint, mInitialTimeBarrier);
    } catch (NuTo::Exception& e)
    {
        e.AddMessage("The simulation failed. \n");
        if (mContinueOnException)
        {
            std::cout << e.ErrorMessage() << "\n but I'll continue.";
        }
        else
        {
            throw;
        }
    }

    std::cout << "min Dist.= " << mParticleHandler->GetAbsoluteMininimalDistance(*mSpecimen)<< std::endl;

    // the spheres are stored with the diameters _after_ the EDMD simulation

}



void NuTo::GeometryConcrete::ExportGmshGeo2D(std::string rGmshFile, double rMeshSize, double rZSlice, double rMinRadius)
{
    if (mParticleHandler == nullptr)
        throw NuTo::Exception("[NuTo::GeometryConcrete::ExportGmsh2D] Run a simulation first!");

    double zS = mSpecimen->GetBoundingBox()(2,0);
    double zE = mSpecimen->GetBoundingBox()(2,1);

    if (zS >= rZSlice || zE <= rZSlice)
        throw NuTo::Exception("[NuTo::GeometryConcrete::ExportGmsh2D] The rZSlice does not intersect the specimen!");

    mParticleHandler->ExportParticlesToGmsh2D(rGmshFile+".geo", *mSpecimen, rMeshSize, rZSlice, rMinRadius);
}

//! @brief ... exports the geometry to a 3D mesh file
//! @param rGmshFile ... path (without .geo) for the gmsh file
//! @param rMeshSize ... mesh size parameter
void NuTo::GeometryConcrete::ExportGmshGeo3D(std::string rGmshFile, double rMeshSize)
{
    if (mParticleHandler == nullptr)
        throw NuTo::Exception("[NuTo::GeometryConcrete::ExportGmsh3D] Run a simulation first!");

    mParticleHandler->ExportParticlesToGmsh3D(rGmshFile+".geo", *mSpecimen, rMeshSize);
}

void NuTo::GeometryConcrete::SetSpecimenBox(double rXs, double rXe, double rYs, double rYe, double rZs, double rZe)
{
    if (mSpecimen)
        delete mSpecimen;

    Eigen::MatrixXd bounds(3,2);
    bounds << rXs, rXe, rYs, rYe, rZs, rZe;
    mSpecimen = new NuTo::Specimen(bounds, NuTo::Specimen::Box);
}

void NuTo::GeometryConcrete::SetSpecimenCylinder(double rXs, double rXe, double rYs, double rYe, double rZs, double rZe)
{
    Eigen::MatrixXd bounds(3,2);
    bounds << rXs, rXe, rYs, rYe, rZs, rZe;
    mSpecimen = new NuTo::Specimen(bounds, NuTo::Specimen::Cylinder);
}

void NuTo::GeometryConcrete::SetGradingCurve(eGradingCurve rGradingCurveEnum, int rNumClasses)
{
    Eigen::MatrixX3d fullGradingCurve(6,3);
    switch (rGradingCurveEnum)
    {
        case A16:
            fullGradingCurve <<
              8,   16, 0.40,
              4,    8, 0.24,
              2,    4, 0.15,
              1,    2, 0.09,
              0.5,  1, 0.04,
              0.25,.5, 0.05;
            break;
        case B16:
            fullGradingCurve <<
              8,   16, 0.24,
              4,    8, 0.20,
              2,    4, 0.14,
              1,    2, 0.10,
              0.5,  1, 0.12,
              0.25,.5, 0.12;
            break;
        case C16:
            fullGradingCurve <<
              8,   16, 0.12,
              4,    8, 0.14,
              2,    4, 0.12,
              1,    2, 0.13,
              0.5,  1, 0.15,
              0.25,.5, 0.16;
            break;
        default:
            throw NuTo::Exception("[NuTo::GeometryConcrete::SetGradingCurve] Desired type currently not implemented.");
            break;
    }

    if (rNumClasses < 1 || rNumClasses > 6)
        throw NuTo::Exception("[NuTo::GeometryConcrete::SetGradingCurve] 1 < rNumClasses <= 6 !");

    mGradingCurve = fullGradingCurve.block(0,0,rNumClasses, 3);
    std::cout << "[NuTo::GeometryConcrete::SetGradingCurve] Used grading curve: \n" << mGradingCurve << std::endl;
}

void NuTo::GeometryConcrete::SetGradingCurve(const Eigen::MatrixXd& rGradingCurve)
{
    mGradingCurve = rGradingCurve;
}



void NuTo::GeometryConcrete::CheckParameters()
{
    if (mSpecimen == nullptr)
        throw NuTo::Exception("[NuTo::GeometryConcrete::CheckParameters] Specimen not defined. Set it first");

    if (mGradingCurve.size() == 0)
        throw NuTo::Exception("[NuTo::GeometryConcrete::CheckParameters] Grading curve not defined. Set it first");

}

Eigen::MatrixXd NuTo::GeometryConcrete::GetParticles(bool rBeforeEDMD)
{
    if (mParticleHandler == nullptr)
        throw NuTo::Exception("[NuTo::GeometryConcrete::GetSpheres] Run a simulation first!");

    return mParticleHandler->GetParticles(rBeforeEDMD);
}

void NuTo::GeometryConcrete::SetParticles(Eigen::MatrixXd rParticles)
{
    mParticleHandler = new NuTo::ParticleHandler(rParticles, mRandomVelocityRange, mRelativeGrowthRate, mAbsoluteGrowthRate);
}


void NuTo::GeometryConcrete::SetAbsoluteGrowthRate(double rAbsoluteGrowthRate)
{
    mAbsoluteGrowthRate = rAbsoluteGrowthRate;
    mRelativeGrowthRate = 0.;
}

void NuTo::GeometryConcrete::SetRelativeGrowthRate(double rRelativeGrowthRate)
{
    mRelativeGrowthRate = rRelativeGrowthRate;
    mAbsoluteGrowthRate = 0.;
}

void NuTo::GeometryConcrete::SetInitialTimeBarrier(double rInitialTimeBarrier)
{
    mInitialTimeBarrier = rInitialTimeBarrier;
}

void NuTo::GeometryConcrete::SetNumEventsMax(double rNumEventsMax)
{
    mNumEventsMax = rNumEventsMax;
}

void NuTo::GeometryConcrete::SetParticleVolumeFraction(double rParticleVolumeFraction)
{
    mParticleVolumeFraction = rParticleVolumeFraction;
}

void NuTo::GeometryConcrete::SetRandomVelocityRange(double rRandomVelocityRange)
{
    mRandomVelocityRange = rRandomVelocityRange;
}

void NuTo::GeometryConcrete::SetSecondsPrint(double rSecondsPrint)
{
    mSecondsPrint = rSecondsPrint;
}

void NuTo::GeometryConcrete::SetSecondsWallTimeMax(double rSecondsWallTimeMax)
{
    mSecondsWallTimeMax = rSecondsWallTimeMax;
}

void NuTo::GeometryConcrete::SetSeed(double rSeed)
{
    mSeed = rSeed;
}

bool NuTo::GeometryConcrete::ContinueOnException() const
{
    return mContinueOnException;
}

void NuTo::GeometryConcrete::SetContinueOnException(bool rContinueOnException)
{
    mContinueOnException = rContinueOnException;
}
