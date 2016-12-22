%module(package="nuto") ModulGeometryConcrete
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/FullMatrix.h"
#include "math/FullVector.h"
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

%apply int& OUTPUT { int& rType };
double NuTo::CollidableParticleSphere::PredictCollision(NuTo::CollidableParticleSphere& rSphere, int& rType);
double NuTo::CollidableParticleSphere::PredictCollision(NuTo::CollidableWallPhysical& rWall, int& rType);


// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
%ignore Exception;
%include "base/ModulNuToBase.i"

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


