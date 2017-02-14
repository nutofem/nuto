# -*- coding: utf-8 -*-

import nuto
import numpy as np
import unittest


class Brick8NTestCase(unittest.TestCase):
    def test_SphereSphere(self):
        sphere1 = nuto.CollidableParticleSphere(np.array([0.,0.,0.]), np.array([3.,0.,0.]), 2., 1., 10)
        # surface position after timeCollision = 2 seconds: x = 0 + 3*2 + 2 + 1*2 = 10

        sphere2 = nuto.CollidableParticleSphere(np.array([19.,0.,0.]), np.array([-2.,0.,0.]), 1., 2., 11)
        # surface position after timeCollision = 2 seconds: x = 19 - 2*2 - 1 - 2*2 = 10

        timeCollision, eventType = sphere1.PredictCollision(sphere2)
        self.assertAlmostEqual(timeCollision, 2., msg="synchronous")

        sphere1.MoveAndGrow(-42.42)
        sphere2.MoveAndGrow(-6.174)
        timeCollision, eventType = sphere1.PredictCollision(sphere2)
        self.assertAlmostEqual(timeCollision, 2., msg="asynchronous")

        event = nuto.Event(timeCollision, sphere1, sphere2, eventType)
        event.PerformCollision()
        timeCollision, eventType = sphere1.PredictCollision(sphere2)
        self.assertAlmostEqual(timeCollision, -1, msg="after collision")

    def test_SphereSphereGrowthOnly(self):
        sphere1 = nuto.CollidableParticleSphere(np.array([0.,0.,0.]), np.array([0.,0.,0.]), 2., 1., 10)
        # surface position after timeCollision = 1 second: x = 0 + 2 + 1*1 = 3

        sphere2 = nuto.CollidableParticleSphere(np.array([6.,0.,0.]), np.array([0.,0.,0.]), 1., 2., 11)
        # surface position after timeCollision = 1 second: x = 6 - 1 - 2*1 = 2

        sphere1.MoveAndGrow(.5)
        sphere2.MoveAndGrow(.2)

        timeCollision, eventType = sphere1.PredictCollision(sphere2)
        self.assertAlmostEqual(timeCollision, 1, msg="CollisionGrowthOnlySphere")

    def test_SphereWall(self):
        sphere1 = nuto.CollidableParticleSphere(np.array([0.,0.,0.]), np.array([3.,0.,0.]), 2., 1., 10)
            # surface position after timeCollision = 2 seconds: x = 0 + 3*2 + 2 + 1*2 = 10
        wall = nuto.CollidableWallPhysical(np.array([10.,0.,0.]), np.array([-1.,0.,0.]), 0)

        timeCollision, eventType = sphere1.PredictCollision(wall)
        self.assertAlmostEqual(timeCollision, 2, msg="synchronous")

        event = nuto.Event(timeCollision, sphere1, wall, eventType)
        event.PerformCollision()
        timeCollision, eventType = sphere1.PredictCollision(wall)
        self.assertAlmostEqual(timeCollision, -1, msg="after collision")

    def test_EventList(self):

        sphere1 = nuto.CollidableParticleSphere(np.array([0., 0.,0.]), np.array([0.,0.,0.]), 0, 0, 1)
        sphere2 = nuto.CollidableParticleSphere(np.array([0., 0.,0.]), np.array([0.,0.,0.]), 0, 0, 1)
        sphere3 = nuto.CollidableParticleSphere(np.array([0., 0.,0.]), np.array([0.,0.,0.]), 0, 0, 1)
        sphere4 = nuto.CollidableParticleSphere(np.array([0., 0.,0.]), np.array([0.,0.,0.]), 0, 0, 1)

        events = nuto.EventListHandler()
        events.AddEvent(0., sphere1, sphere2, 0)
        events.AddEvent(0., sphere2, sphere1, 0)
        events.AddEvent(0., sphere1, sphere2, 0)
        self.assertEqual(events.GetEventListSize(), 1, msg="simultaneous, same events")
        events.Clear()

        events.AddEvent(0., sphere2, sphere1, 0)
        events.AddEvent(0., sphere3, sphere1, 0)
        events.AddEvent(0., sphere4, sphere1, 0)
        self.assertEqual(events.GetEventListSize(), 3, msg="simultaneous, differ first")
        events.Clear()

        events.AddEvent(0., sphere1, sphere2, 0)
        events.AddEvent(0., sphere1, sphere3, 0)
        events.AddEvent(0., sphere1, sphere4, 0)
        self.assertEqual(events.GetEventListSize(), 3, msg="simultaneous, differ first")

    def test_Box(self):
        numParticles = 1000
        bBoxLength = 40.
        bBox = np.array([[-bBoxLength/2.,bBoxLength/2.],
                         [-bBoxLength/2.,bBoxLength/2.],
                         [-bBoxLength/2.,bBoxLength/2.]])

        specimen = nuto.Specimen(bBox, 0)

        spheres = nuto.ParticleHandler(numParticles, bBox, 1., 0.1)

        subBoxes = nuto.SubBoxHandler(spheres, specimen, 10)

        collisions = nuto.CollisionHandler(spheres, subBoxes, "")

        collisions.Simulate(
           100000000,
           8.,
           30.,
           1.,
           10.)

    def test_Cylinder(self):

        numParticles = 1000
        bBoxLength = 40. * .7
        bBoxCylPos = np.array([[-bBoxLength/2.,bBoxLength/2.],
                               [-bBoxLength/2.,bBoxLength/2.],
                               [-bBoxLength/2.,bBoxLength/2.]])

        specimen = nuto.Specimen(bBoxCylPos, 2)

        spheres = nuto.ParticleHandler(numParticles, bBoxCylPos, 1., 0.1)

        subBoxes = nuto.SubBoxHandler(spheres, specimen, 10)

        collisions = nuto.CollisionHandler(spheres, subBoxes, "")

        collisions.Simulate(
           100000000,
           8.,
           30,
           1.,
           10.)