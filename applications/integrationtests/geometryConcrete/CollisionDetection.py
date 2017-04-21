#!/usr/bin/env python3
import unittest
import numpy as np
import nuto


class Brick8NTestCase(unittest.TestCase):
    def test_SphereSphere(self):
        # surface position after timeCollision = 2 seconds:
        # x = 0 + 3*2 + 2 + 1*2 = 10
        sphere1 = nuto.CollidableParticleSphere(np.r_[0., 0., 0.],
                                                np.r_[3., 0., 0.], 2., 1., 10)

        # surface position after timeCollision = 2 seconds:
        # x = 19 - 2*2 - 1 - 2*2 = 10
        sphere2 = nuto.CollidableParticleSphere(np.r_[19., 0., 0.],
                                                np.r_[-2., 0., 0.], 1., 2., 11)

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
        # surface position after timeCollision = 1 second:
        # x = 0 + 2 + 1*1 = 3
        sphere1 = nuto.CollidableParticleSphere(np.r_[0., 0., 0.],
                                                np.r_[0., 0., 0.], 2., 1., 10)

        # surface position after timeCollision = 1 second:
        # x = 6 - 1 - 2*1 = 2
        sphere2 = nuto.CollidableParticleSphere(np.r_[6., 0., 0.],
                                                np.r_[0., 0., 0.], 1., 2., 11)

        sphere1.MoveAndGrow(.5)
        sphere2.MoveAndGrow(.2)
        timeCollision, eventType = sphere1.PredictCollision(sphere2)
        self.assertAlmostEqual(timeCollision, 1, msg="CollisionGrowthOnlySphere")

    def test_SphereWall(self):
        # surface position after timeCollision = 2 seconds:
        # x = 0 + 3*2 + 2 + 1*2 = 10
        sphere1 = nuto.CollidableParticleSphere(np.r_[0., 0., 0.],
                                                np.r_[3., 0., 0.], 2., 1., 10)
        wall = nuto.CollidableWallPhysical(np.r_[10., 0., 0.],
                                           np.r_[-1., 0., 0.], 0)

        timeCollision, eventType = sphere1.PredictCollision(wall)
        self.assertAlmostEqual(timeCollision, 2, msg="synchronous")

        event = nuto.Event(timeCollision, sphere1, wall, eventType)
        event.PerformCollision()
        timeCollision, eventType = sphere1.PredictCollision(wall)
        self.assertAlmostEqual(timeCollision, -1, msg="after collision")

    def test_EventList(self):
        origin = np.r_[0., 0., 0.]
        sphere1 = nuto.CollidableParticleSphere(origin, origin, 0, 0, 1)
        sphere2 = nuto.CollidableParticleSphere(origin, origin, 0, 0, 1)
        sphere3 = nuto.CollidableParticleSphere(origin, origin, 0, 0, 1)
        sphere4 = nuto.CollidableParticleSphere(origin, origin, 0, 0, 1)

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
        bBox = np.array([[-bBoxLength/2., bBoxLength/2.],
                         [-bBoxLength/2., bBoxLength/2.],
                         [-bBoxLength/2., bBoxLength/2.]])

        specimen = nuto.Specimen(bBox, 0)
        spheres = nuto.ParticleHandler(numParticles, bBox, 1., 0.1)
        subBoxes = nuto.SubBoxHandler(spheres, specimen, 10)
        collisions = nuto.CollisionHandler(spheres, subBoxes, "")

        collisions.Simulate(100000000, 8., 30., 1., 10.)

    def test_Cylinder(self):
        numParticles = 1000
        bBoxLength = 40.0
        bBoxCylPos = np.array([[-bBoxLength/2., bBoxLength/2.],
                               [-bBoxLength/2., bBoxLength/2.],
                               [-bBoxLength/2., bBoxLength/2.]])

        specimen = nuto.Specimen(bBoxCylPos, 2)
        spheres = nuto.ParticleHandler(numParticles, 0.7*bBoxCylPos, 1., 0.1)
        subBoxes = nuto.SubBoxHandler(spheres, specimen, 10)
        collisions = nuto.CollisionHandler(spheres, subBoxes, "")

        collisions.Simulate(100000000, 8., 30, 1., 10.)


if __name__ == "__main__":
    unittest.main()
