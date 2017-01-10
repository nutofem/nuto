%module(package="nuto") ModulGeometryConcrete
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "geometryConcrete/Specimen.h"
#include "geometryConcrete/collision/collidables/CollidableBase.h"
#include "geometryConcrete/collision/collidables/CollidableParticleBase.h"
#include "geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "geometryConcrete/collision/collidables/CollidableWallPhysical.h"
#include "geometryConcrete/collision/Event.h"
#include "geometryConcrete/collision/SubBox.h"
#include "geometryConcrete/collision/handler/SubBoxHandler.h"
#include "geometryConcrete/collision/handler/EventListHandler.h"
#include "geometryConcrete/collision/handler/ParticleHandler.h"
#include "geometryConcrete/collision/handler/CollisionHandler.h"
#include "geometryConcrete/takeAndPlace/ParticleCreator.h"
#include "geometryConcrete/InputReader.h"
%}


%include "math/NuToMath.i" // defines typenames for std::vector and Eigen::Matrix

%apply int& OUTPUT { int& rType };
double NuTo::CollidableParticleSphere::PredictCollision(NuTo::CollidableParticleSphere& rSphere, int& rType);
double NuTo::CollidableParticleSphere::PredictCollision(NuTo::CollidableWallPhysical& rWall, int& rType);



%include "geometryConcrete/Specimen.h"
%include "geometryConcrete/collision/Event.h"
%include "geometryConcrete/collision/SubBox.h"
%include "geometryConcrete/collision/collidables/CollidableBase.h"
%include "geometryConcrete/collision/collidables/CollidableParticleBase.h"
%include "geometryConcrete/collision/collidables/CollidableParticleSphere.h"
%include "geometryConcrete/collision/collidables/CollidableWallBase.h"
%include "geometryConcrete/collision/collidables/CollidableWallPhysical.h"
%include "geometryConcrete/collision/handler/EventListHandler.h"
%include "geometryConcrete/collision/handler/ParticleHandler.h"
%include "geometryConcrete/collision/handler/SubBoxHandler.h"
%include "geometryConcrete/collision/handler/CollisionHandler.h"
%include "geometryConcrete/takeAndPlace/ParticleCreator.h"
%include "geometryConcrete/InputReader.h"


