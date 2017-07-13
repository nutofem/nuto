#include "BoostUnitTest.h"
#include "geometryConcrete/GeometryConcrete.h"
#include "base/Exception.h"

NuTo::GeometryConcrete DefineGeometry(int numClasses, double l = 16)
{
    NuTo::GeometryConcrete geometry;
    geometry.SetSeed(1337);
    geometry.SetSpecimenBox(0, l, 0, l, 0, l);
    geometry.SetGradingCurve(NuTo::GeometryConcrete::B16, numClasses);
    geometry.SetInitialTimeBarrier(10);
    return geometry;
}

BOOST_AUTO_TEST_CASE(API_2D_ParticleDistance)
{
    auto geometry = DefineGeometry(2, 16);
    geometry.SetAbsoluteGrowthRate(0.1);
    geometry.SetParticleVolumeFraction(0.4);
    geometry.MaximizeParticleDistance(0.75);
    geometry.ExportGmshGeo2D("./GeometryConcreteAPI2D", 0.75, 8., 0.75);
}

BOOST_AUTO_TEST_CASE(API_2D_ParticleDistance_Continue)
{
    auto geometry = DefineGeometry(1, 30);
    geometry.SetAbsoluteGrowthRate(0.1);
    geometry.SetParticleVolumeFraction(0.4);
    geometry.SetNumEventsMax(1);
    BOOST_CHECK_THROW(geometry.MaximizeParticleDistance(10), NuTo::Exception);
    geometry.SetContinueOnException(true);
    BOOST_CHECK_NO_THROW(geometry.MaximizeParticleDistance(10));
}

BOOST_AUTO_TEST_CASE(API_3D_VolumeFraction)
{
    auto geometry = DefineGeometry(2, 20);
    geometry.SetParticleVolumeFraction(0.8);
    geometry.SetRelativeGrowthRate(0.1);
    geometry.SetInitialTimeBarrier(2);
    geometry.MaximizeParticleVolumeFraction(0.05);
    geometry.ExportGmshGeo3D("./GeometryConcreteAPI2D", 0.75);
}

BOOST_AUTO_TEST_CASE(API_3D_VolumeFraction_Continue)
{
    auto geometry = DefineGeometry(2, 20);
    geometry.SetParticleVolumeFraction(0.8);
    geometry.SetRelativeGrowthRate(0.1);
    geometry.SetInitialTimeBarrier(2);
    geometry.SetNumEventsMax(1);
    BOOST_CHECK_THROW(geometry.MaximizeParticleVolumeFraction(0.2), NuTo::Exception);
    geometry.SetContinueOnException(true);
    BOOST_CHECK_NO_THROW(geometry.MaximizeParticleVolumeFraction(0.2));
}
