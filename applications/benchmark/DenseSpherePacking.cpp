#include <mechanics/MechanicsEnums.h>
#include "Benchmark.h"
#include "geometryConcrete/collision/handler/CollisionHandler.h"
#include "geometryConcrete/collision/handler/ParticleHandler.h"
#include "geometryConcrete/collision/handler/SubBoxHandler.h"
#include "geometryConcrete/Specimen.h"

BENCHMARK(DenseSpherePacking, Run, runner)
{
    Eigen::MatrixXd boundingBox(3,2);
    boundingBox << 0, 1, 0, 1, 0, 1;


    NuTo::Specimen s(boundingBox, NuTo::Specimen::Box);

    constexpr int numParticles = 1000;
    constexpr double randomVelocityRange = 0.001;
    constexpr double growthRate = 0.01;


    constexpr int numEvents = 10000;
    constexpr double timeMax = 5;
    constexpr double wallTimeMax = 10.;
    constexpr double timePrintout = 1e3;
    constexpr double initialTimeBarrier = 0.01;

    while(runner.KeepRunningIterations(10))
    {
        NuTo::ParticleHandler particleHandler(numParticles, s.GetBoundingBox(), randomVelocityRange, growthRate);
        NuTo::SubBoxHandler subBoxHandler(particleHandler, s, 10);
        NuTo::CollisionHandler collisionHandler(particleHandler, subBoxHandler, "");
        collisionHandler.Simulate(numEvents, timeMax, wallTimeMax, timePrintout, initialTimeBarrier);
    }
}