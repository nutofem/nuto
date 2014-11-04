# -*- coding: utf-8 -*-

import math
# load nuto package
import nuto
import sys


error = False


# ======================================================
# ==       Sphere-Sphere collision detection          ==
# ======================================================

sphere1 = nuto.CollidableParticleSphere(nuto.DoubleFullVector((0.,0.,0.)), nuto.DoubleFullVector((3.,0.,0.)), 2., 1., 10)
    # surface position after timeCollision = 2 seconds: x = 0 + 3*2 + 2 + 1*2 = 10

sphere2 = nuto.CollidableParticleSphere(nuto.DoubleFullVector((19.,0.,0.)), nuto.DoubleFullVector((-2.,0.,0.)), 1., 2., 11)
    # surface position after timeCollision = 2 seconds: x = 19 - 2*2 - 1 - 2*2 = 10

    
print "[CollisionFrontalSphere] synchronous"
timeCollision, eventType = sphere1.PredictCollision(sphere2)
if abs(timeCollision - 2.) > 1e-14 :
    print "Collision Time = ", timeCollision, "  E R R O R"
    error = True
  
  
print "[CollisionFrontalSphere] asynchronous"
sphere1.MoveAndGrow(-42.42)
sphere2.MoveAndGrow(-6.174)
timeCollision, eventType = sphere1.PredictCollision(sphere2)
if abs(timeCollision - 2.) > 1e-14 :
    print "Collision Time = ", timeCollision, "  E R R O R"
    error = True

print "[CollisionFrontalSphere] after collision"
event = nuto.Event(timeCollision, sphere1, sphere2, eventType)
event.PerformCollision()
timeCollision, eventType = sphere1.PredictCollision(sphere2)
if timeCollision != -1 > 1e-14 :
    print "Collision Time = ", timeCollision, "  E R R O R"
    error = True

# ======================================================
# ==  Sphere-Sphere collision detection, growth only  ==
# ======================================================  

sphere1 = nuto.CollidableParticleSphere(nuto.DoubleFullVector((0.,0.,0.)), nuto.DoubleFullVector((0.,0.,0.)), 2., 1., 10)
    # surface position after timeCollision = 1 second: x = 0 + 2 + 1*1 = 3
sphere2 = nuto.CollidableParticleSphere(nuto.DoubleFullVector((6.,0.,0.)), nuto.DoubleFullVector((0.,0.,0.)), 1., 2., 11)
    # surface position after timeCollision = 1 second: x = 6 - 1 - 2*1 = 2
sphere1.MoveAndGrow(.5)
sphere2.MoveAndGrow(.2)
    
print "[CollisionGrowthOnlySphere]"
timeCollision, eventType = sphere1.PredictCollision(sphere2)
if abs(timeCollision - 1.) > 1e-14 :
    print "Collision Time = ", timeCollision, "  E R R O R"
    error = True
    
    
    
# ======================================================
# ==        Wall-Sphere collision detection           ==
# ======================================================

sphere1 = nuto.CollidableParticleSphere(nuto.DoubleFullVector((0.,0.,0.)), nuto.DoubleFullVector((3.,0.,0.)), 2., 1., 10)
    # surface position after timeCollision = 2 seconds: x = 0 + 3*2 + 2 + 1*2 = 10

wall = nuto.CollidableWallPhysical(nuto.DoubleFullVector((10.,0.,0.)), nuto.DoubleFullVector((-1.,0.,0.)), 0)
    
print "[CollisionFrontalWall] synchronous"
timeCollision, eventType = sphere1.PredictCollision(wall)
if abs(timeCollision - 2.) > 1e-14 :
    print "Collision Time = ", timeCollision, "  E R R O R"
    error = True

print "[CollisionFrontalWall] after collision"
event = nuto.Event(timeCollision, sphere1, wall, eventType)
event.PerformCollision()
timeCollision, eventType = sphere1.PredictCollision(wall)
if timeCollision != -1 > 1e-14 :
    print "Collision Time = ", timeCollision, "  E R R O R"
    error = True

# ======================================================
# ==       simultaneous collision detection           ==
# ======================================================

# build sphere matrix
#rawSpheres = nuto.DoubleFullMatrix(3,4, [0.,0.,0.,4,           10.,0.,0.,4.,            5.,math.sqrt(75.),0.,4.])
rawSpheres = nuto.DoubleFullMatrix(3,4, [0.,10.,5.,  0.,0.,math.sqrt(75.),   0.,0.,0., 5.,5.,5.])


# build boundary matrix
boundary = nuto.DoubleFullMatrix(3,2, [-1e10, -1e10, -1e10, 1e10, 1e10, 1e10])
specimen = nuto.Specimen(boundary, 0);

# build sub box division vector
subBoxDivs = nuto.IntFullVector([1,1,1])

# particle handler
spheres = nuto.ParticleHandler(rawSpheres,0.,0.,1.)

# sub box handler
subBoxes = nuto.SubBoxHandler(spheres, specimen, subBoxDivs)
subBoxes.Build()

events = nuto.EventListHandler()
dummy = events.SetTimeBarrier(1000.,subBoxes)

#.........
 
 
# ======================================================
# ==       EventList behaviour                        ==
# ====================================================== 


sphere3 = nuto.CollidableParticleSphere(nuto.DoubleFullVector((0., 0.,0.)), nuto.DoubleFullVector((0.,0.,0.)), 0, 0, 1)
sphere4 = nuto.CollidableParticleSphere(nuto.DoubleFullVector((0., 0.,0.)), nuto.DoubleFullVector((0.,0.,0.)), 0, 0, 1)

events = nuto.EventListHandler()
print ""
print "[EventList behaviour] simultaneous, same events"
events.AddEvent(0., sphere1, sphere2, 0)
events.AddEvent(0., sphere2, sphere1, 0)
events.AddEvent(0., sphere1, sphere2, 0)
numEvents = events.GetEventListSize();
if (events.GetEventListSize == 1) :
  print " E R R O R", events.GetEventListSize()
  error = True
events.Clear()    
 
print "[EventList behaviour] simultaneous, differ first"
events.AddEvent(0., sphere2, sphere1, 0)
events.AddEvent(0., sphere3, sphere1, 0)
events.AddEvent(0., sphere4, sphere1, 0)
if events.GetEventListSize == 3 :
  print " E R R O R"
  error = True
events.Clear()

print "[EventList behaviour] simultaneous, differ second"
events.AddEvent(0., sphere1, sphere2, 0)
events.AddEvent(0., sphere1, sphere3, 0)
events.AddEvent(0., sphere1, sphere4, 0)
if events.GetEventListSize == 3 :
  print " E R R O R"
  error = True
events.Clear()
 
events.Clear()
 
 
# ======================================================
# ==      Complex test of a box                       ==
# ====================================================== 
 
numParticles = 1000;
bBoxLength = 40.
bBox = nuto.DoubleFullMatrix(3,2,[-bBoxLength/2., -bBoxLength/2., -bBoxLength/2., bBoxLength/2., bBoxLength/2., bBoxLength/2.])

specimen = nuto.Specimen(bBox, 0)

spheres = nuto.ParticleHandler(numParticles, bBox, 1., 0.1)

subBoxes = nuto.SubBoxHandler(spheres, specimen, 10)
subBoxes.Build(4)

collisions = nuto.CollisionHandler(spheres, subBoxes, "")

collisions.Simulate(
   100000000,
   8.,
   30.,
   1.,
   10.) 
spheres = None
subBoxes = None
collisions = None

   
# ======================================================
# ==      Complex test of a cylinder                  ==
# ====================================================== 


bBoxLength = 40. * .7
bBoxCylPos = nuto.DoubleFullMatrix(3,2,[-bBoxLength/2., -bBoxLength/2., -bBoxLength/2., bBoxLength/2., bBoxLength/2., bBoxLength/2.])

specimen = nuto.Specimen(bBox, 2)

spheres = nuto.ParticleHandler(numParticles, bBoxCylPos, 1., 0.1)

subBoxes = nuto.SubBoxHandler(spheres, specimen, 10)
subBoxes.Build(4)

collisions = nuto.CollisionHandler(spheres, subBoxes, "")

collisions.Simulate(
   100000000,
   8.,
   30,
   1.,
   10.) 
 

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
