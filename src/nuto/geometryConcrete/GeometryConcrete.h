/*
 * GeometryConcrete.h
 *
 *  Created on: 4 Sep 2015
 *      Author: ttitsche
 */

#ifndef GEOMETRYCONCRETE_GEOMETRYCONCRETE_H_
#define GEOMETRYCONCRETE_GEOMETRYCONCRETE_H_

#include "nuto/math/FullMatrix_Def.h"

namespace NuTo
{

class Specimen;
class ParticleHandler;

class GeometryConcrete
{
public:

    enum eGradingCurve
    {
        A16,
        B16,
        C16
    };

    GeometryConcrete();

    virtual ~GeometryConcrete();


    //! @brief runs a RSA simulation with the desired phi
    //! afterwards, runs a EDMD simulation to maximize the particle distance
    //! @param rParticleDistance ... desired particle distance
    void MaximizeParticleDistance(double rParticleDistance);


    //! @brief ... runs a RSA simulation with smaller spheres
    //! afterwards, runs a EDMD simulation to reach the original sphere size
    //! @param rShrinkage ... sphere distribution gets multiplied by (1-rShrinkage)
    void MaximizeParticleVolumeFraction(double rShrinkage);


    //! @brief ... exports the geometry to a 2D mesh file
    //! @param rGmshFile ... path (without .geo) for the gmsh file
    //! @param rMeshSize ... mesh size parameter
    //! @param rZSlice ... 3D --> 2D at z=rZSlice
    void ExportGmshGeo2D(std::string rGmshFile, double rMeshSize, double rZSlice);


    //! @brief ... exports the geometry to a 3D mesh file
    //! @param rGmshFile ... path (without .geo) for the gmsh file
    //! @param rMeshSize ... mesh size parameter
    void ExportGmshGeo3D(std::string rGmshFile, double rMeshSize);


    void SetSpecimenBox(double rXs, double rXe, double rYs, double rYe, double rZs, double rZe);
    void SetSpecimenCylinder(double rXs, double rXe, double rYs, double rYe, double rZs, double rZe);

    //! @brief ... sets the grading curve - DIN
    //! @param rGradingCurveEnum ... enum
    //! @param rNumClasses ... number of classes , beginning from the biggest particles, to consider
    void SetGradingCurve(eGradingCurve rGradingCurveEnum, int rNumClasses);

    //! @brief ... sets the grading curve - self defined
    //! @param rGradingCurve ... user defined grading curve
    void SetGradingCurve(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rGradingCurve);


    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> GetParticles(bool rBeforeEDMD);

    void SetParticles(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>);

    //! @brief sets the absolute grow rate, removes the relative growth rate
    void SetAbsoluteGrowthRate(double rAbsoluteGrowthRate);

    //! @brief sets the relative grow rate, removes the absolute growth rate
    void SetRelativeGrowthRate(double rRelativeGrowthRate);

    void SetGmshBinary(const std::string& rGmshBinary);

    void SetInitialTimeBarrier(double rInitialTimeBarrier);

    void SetNumEventsMax(double rNumEventsMax);

    void SetParticleVolumeFraction(double rParticleVolumeFraction);

    void SetRandomVelocityRange(double rRandomVelocityRange);

    void SetSecondsPrint(double rSecondsPrint);

    void SetSecondsWallTimeMax(double rSecondsWallTimeMax);

    void SetSeed(double rSeed);

private:

    void CheckParameters();

    Specimen* mSpecimen = nullptr;
    ParticleHandler* mParticleHandler = nullptr;
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mGradingCurve;

    double mSeed = 42;
    double mParticleVolumeFraction = 0.5;

    double mNumEventsMax = 2.e6;
    double mSecondsWallTimeMax = 60.;
    double mSecondsPrint = 5.;
    double mInitialTimeBarrier = 10;

    double mRandomVelocityRange = 1.;
    double mRelativeGrowthRate = 0.0;
    double mAbsoluteGrowthRate = 0.0;

};

} /* namespace NuTo */

#endif /* GEOMETRYCONCRETE_GEOMETRYCONCRETE_H_ */
