/*
 * CollidableWallCylinder.cpp
 *
 *  Created on: 13 Feb 2014
 *      Author: ttitsche
 */

#include <eigen3/Eigen/Dense> // for cross product

#include "base/Exception.h"
#include "geometryConcrete/collision/Event.h"
#include "geometryConcrete/collision/collidables/CollidableWallCylinder.h"
#include "geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "visualize/UnstructuredGrid.h"

NuTo::CollidableWallCylinder::CollidableWallCylinder(Eigen::Vector3d rPosition, Eigen::Vector3d rDirection,
                                                     double rRadius, double rHeigth, int rIndex)
    : CollidableWallBase(rPosition, rDirection, rIndex)
    , mRadius(rRadius)
    , mHeigth(rHeigth)
{
}

void NuTo::CollidableWallCylinder::PerformCollision(CollidableParticleSphere& rSphere)
{
    // get perpendicular from sphere centre to direction vector
    Eigen::Vector3d dP = mPosition - rSphere.mPosition;
    Eigen::Vector3d n = dP - dP.dot(mDirection) * mDirection;
    if (n.dot(n) == 0)
        throw NuTo::Exception("[NuTo::CollidableWallCylinder::PerformCollision] n = 0");

    n.normalize();

    // perform plane wall collision:
    double velocityNormal = rSphere.mVelocity.dot(n);

    if (velocityNormal < 0.)
        rSphere.mVelocity += -2 * velocityNormal * n;

    rSphere.mVelocity += (rSphere.mGrowthRate) * n;
}

double NuTo::CollidableWallCylinder::PredictCollision(CollidableParticleSphere& rSphere, int& rType)
{
    rType = Event::EventType::WallCollision;

    Eigen::Vector3d dPos = rSphere.mPosition - mPosition;

    Eigen::Vector3d dP = dPos - dPos.dot(mDirection) * mDirection;
    Eigen::Vector3d dV = rSphere.mVelocity - rSphere.mVelocity.dot(mDirection) * mDirection;
    long double dR = mRadius - rSphere.mRadius;
    long double dG = -rSphere.mGrowthRate;

    long double a = dV.dot(dV) - dG * dG;
    long double c = dP.dot(dP) - dR * dR;
    long double b = 2 * (dV.dot(dP) - dR * dG);


    if (c > 2.0e-10 * rSphere.mRadius)
        throw NuTo::Exception("[NuTo::CollidableSphere::CreateNewEvent] Sphere overlapping with Cylinder!");

    long double discriminant = b * b - 4 * a * c;
    if (discriminant > -1e-11)
    {
        long double timeCollision = 0;
        if (discriminant < 0.)
            discriminant = 0;

        if (b > 0.)
        {

            timeCollision = 2 * c / (-b - std::sqrt(discriminant));
            if (timeCollision >= 0.)
                return rSphere.mTimeOfLastUpdate + timeCollision;
        }
        else if (a > 0. && b <= 0.)
        {
            timeCollision = (-b + std::sqrt(discriminant)) / (2 * a);
            return rSphere.mTimeOfLastUpdate + timeCollision;
        }
    }
    return Event::EVENTNULL;
}

void NuTo::CollidableWallCylinder::VisualizationStatic(Visualize::UnstructuredGrid& rVisualizer) const
{
    // =========================================
    // create a circle out of several line cells
    // =========================================
    const int nPoints = 40; // resolution of the circle;

    Eigen::Vector3d n = mDirection; // rename
    Eigen::Vector3d nonCoLin({1, 0, 0});
    if (nonCoLin == n)
        nonCoLin << 0, 1, 0;

    Eigen::Vector3d sVec = n.cross(nonCoLin); // random vector perpendicular to mDirection = n
    sVec.normalize();

    Eigen::Matrix<double, nPoints + 1, 3> pointCircle; // store circle coordinates
    double dAlpha = 2. * M_PI / nPoints;
    for (int i = 0; i < nPoints + 1; i++)
    {
        // calculate circle points and add them to points
        // point at nPoints+1 closes the circle
        double alpha = i * dAlpha;
        Eigen::Vector3d rotVec =
                n * (n.dot(sVec)) + std::cos(alpha) * n.cross(sVec).cross(n) + std::sin(alpha) * n.cross(sVec);
        pointCircle.row(i) = (mRadius * rotVec).transpose();
    }

    // =========================================
    //   create and add top and bottom points
    // =========================================
    Eigen::Vector3d pointTop = mPosition + n * mHeigth / 2.;
    Eigen::Vector3d pointBot = mPosition - n * mHeigth / 2.;
    unsigned int indexTop = rVisualizer.AddPoint(pointTop);
    unsigned int indexBot = rVisualizer.AddPoint(pointBot);

    unsigned int pointCircleIndexTop[nPoints + 1];
    unsigned int pointCircleIndexBot[nPoints + 1];
    for (int i = 0; i < nPoints + 1; i++)
    {
        Eigen::Vector3d pointCircleTop = pointCircle.row(i).transpose() + pointTop;
        Eigen::Vector3d pointCircleBot = pointCircle.row(i).transpose() + pointBot;
        pointCircleIndexTop[i] = rVisualizer.AddPoint(pointCircleTop);
        pointCircleIndexBot[i] = rVisualizer.AddPoint(pointCircleBot);
    }

    // =========================================
    //        create top and bottom planes
    // =========================================
    for (int i = 0; i < nPoints; i++)
    {
        std::vector<int> triangleIndexTop(3);
        std::vector<int> triangleIndexBot(3);

        triangleIndexTop[0] = indexTop;
        triangleIndexTop[1] = pointCircleIndexTop[i];
        triangleIndexTop[2] = pointCircleIndexTop[i + 1];

        triangleIndexBot[0] = indexBot;
        triangleIndexBot[1] = pointCircleIndexBot[i];
        triangleIndexBot[2] = pointCircleIndexBot[i + 1];

        unsigned int insertIndexTop = rVisualizer.AddCell(triangleIndexTop, eCellTypes::TRIANGLE);
        unsigned int insertIndexBot = rVisualizer.AddCell(triangleIndexBot, eCellTypes::TRIANGLE);

        Eigen::Vector3d dirTop = -n;
        Eigen::Vector3d dirBot = n;
        rVisualizer.SetCellData(insertIndexTop, "Direction", dirTop);
        rVisualizer.SetCellData(insertIndexBot, "Direction", dirBot);
    }
    // =========================================
    //            create side surface
    // =========================================
    for (int i = 0; i < nPoints; i++)
    {
        std::vector<int> quadIndex(4);
        quadIndex[0] = pointCircleIndexBot[i];
        quadIndex[1] = pointCircleIndexBot[i + 1];
        quadIndex[2] = pointCircleIndexTop[i + 1];
        quadIndex[3] = pointCircleIndexTop[i];

        unsigned int insertIndex = rVisualizer.AddCell(quadIndex, eCellTypes::QUAD);
        Eigen::Vector3d lineVector = pointCircle.row(i).transpose() - pointCircle.row(i + 1).transpose();
        Eigen::Vector3d dir = lineVector.cross(n);
        dir.normalize();
        rVisualizer.SetCellData(insertIndex, "Direction", dir);
    }
}


bool NuTo::CollidableWallCylinder::IsInside(const CollidableParticleSphere& rSphere) const
{
    Eigen::Vector3d dPos = rSphere.mPosition - mPosition;
    Eigen::Vector3d dP = dPos - dPos.dot(mDirection) * mDirection;

    return dP.dot(dP) < (mRadius - rSphere.mRadius) * (mRadius - rSphere.mRadius);
}


void NuTo::CollidableWallCylinder::Print(std::ostream&) const
{
    // TODO: NuTo::CollidableWallCylinder::Print
}
