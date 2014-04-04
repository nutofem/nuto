%module(package="nuto") ModulGeometryConcrete
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleBase.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallPhysical.h"
#include "nuto/geometryConcrete/collision/Event.h"
#include "nuto/geometryConcrete/collision/SubBox.h"
#include "nuto/geometryConcrete/collision/handler/SubBoxHandler.h"
#include "nuto/geometryConcrete/collision/handler/EventListHandler.h"
#include "nuto/geometryConcrete/collision/handler/ParticleHandler.h"
#include "nuto/geometryConcrete/collision/handler/CollisionHandler.h"
#include "nuto/geometryConcrete/takeAndPlace/ParticleCreator.h"
#include "nuto/geometryConcrete/InputReader.h"
%}

%apply int& OUTPUT { int& rType };
double NuTo::CollidableParticleSphere::PredictCollision(NuTo::CollidableParticleSphere& rSphere, int& rType);
double NuTo::CollidableParticleSphere::PredictCollision(NuTo::CollidableWallPhysical& rWall, int& rType);


// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
%ignore Exception;
%include "nuto/base/ModulNuToBase.i"


%include "nuto/geometryConcrete/collision/Event.h"
%include "nuto/geometryConcrete/collision/SubBox.h"
%include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"
%include "nuto/geometryConcrete/collision/collidables/CollidableParticleBase.h"
%include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"
%include "nuto/geometryConcrete/collision/collidables/CollidableWallBase.h"
%include "nuto/geometryConcrete/collision/collidables/CollidableWallPhysical.h"
%include "nuto/geometryConcrete/collision/handler/EventListHandler.h"
%include "nuto/geometryConcrete/collision/handler/ParticleHandler.h"
%include "nuto/geometryConcrete/collision/handler/SubBoxHandler.h"
%include "nuto/geometryConcrete/collision/handler/CollisionHandler.h"
%include "nuto/geometryConcrete/takeAndPlace/ParticleCreator.h"
%include "nuto/geometryConcrete/InputReader.h"


