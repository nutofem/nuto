/*
 * SubBoxBuilder.cpp
 *
 *  Created on: 4 Mar 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/WallTime.h"
#include "nuto/geometryConcrete/collision/handler/SubBoxHandler.h"
#include "nuto/geometryConcrete/collision/handler/ParticleHandler.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallBase.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallPhysical.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallVirtual.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallCylinder.h"

#include "nuto/visualize/VisualizeUnstructuredGrid.h"

NuTo::SubBoxHandler::SubBoxHandler(
		ParticleHandler& rSpheres,
		Specimen& rSpecimen,
		const FullVector<int, Eigen::Dynamic> rDivisions)
		:
				mSpecimen(rSpecimen),
				mDivisions(rDivisions)

{
	if (rDivisions.GetNumRows() != 3)
		throw Exception("[NuTo::SubBoxHandler::SubBoxHandler] FullVector<int, 3> rDivisions!");

	mSpheres = &rSpheres;
	mSubBoxes.clear();
}

NuTo::SubBoxHandler::SubBoxHandler(
		ParticleHandler& rSpheres,
		Specimen& rSpecimen,
		const int rParticlesPerSubBox)
		:
				mSpecimen(rSpecimen),
				mDivisions(rSpheres.GetSubBoxDivisions(mSpecimen, rParticlesPerSubBox))
{
	mSpheres = &rSpheres;
	mSubBoxes.clear();

	int numDivisions = mDivisions[0] * mDivisions[1] * mDivisions[2];

	std::cout << "[NuTo::SubBoxHandler::SubBoxHandler] Build " << numDivisions << " sub boxes." << std::endl;

}

NuTo::SubBoxHandler::~SubBoxHandler()
{
	for (auto box : mSubBoxes)
		delete box;
}

unsigned int NuTo::SubBoxHandler::GetBoxIndex(int rX, int rY, int rZ)
{
	// get rX, rY, rZ in 0 .. rDivision - range

	rX = std::min(rX, static_cast<int>(mDivisions[0] - 1));
	rX = std::max(rX, 0);

	rY = std::min(rY, static_cast<int>(mDivisions[1] - 1));
	rY = std::max(rY, 0);

	rZ = std::min(rZ, static_cast<int>(mDivisions[2] - 1));
	rZ = std::max(rZ, 0);

	return mDivisions[1] * mDivisions[2] * rX + mDivisions[2] * rY + rZ;
}

void NuTo::SubBoxHandler::Build(const int rNumThreads)
{
	switch (mSpecimen.GetTypeOfSpecimen())
	{
	case 0:
		BuildBox(rNumThreads);
		break;
	case 2:
		BuildCylinder(rNumThreads);
		break;
	default:
		throw Exception("[NuTo::SubBoxHandler::Build] Type not specimen not implemented yet.");
		break;
	}
}

void NuTo::SubBoxHandler::BuildBox(const int rNumThreads)
{
	std::cout << "[BuildBox] Building boxes..." << std::endl;
	BuildSubBoxes(rNumThreads);
	AddSpheresToBoxes();
}

NuTo::FullVector<double, Eigen::Dynamic> NuTo::SubBoxHandler::GetSubBoxLength()
{
	FullVector<double, 3> subBoxLength = mSpecimen.GetLength();
	subBoxLength[0] /= mDivisions[0];
	subBoxLength[1] /= mDivisions[1];
	subBoxLength[2] /= mDivisions[2];
	return subBoxLength;
}

void NuTo::SubBoxHandler::BuildSubBoxes(const int rNumThreads)
{

	mSubBoxes.clear();

	// equal length sub-boxes
	FullVector<double, 3> subBoxLength = GetSubBoxLength();
	int index = 0;

	unsigned int nDivisions = mDivisions[0] * mDivisions[1] * mDivisions[2];

	// create cell grid boxes
	for (unsigned int i = 0; i < nDivisions; ++i)
		mSubBoxes.push_back(new SubBox(i, rNumThreads));

	// TODO make this ugly method readable, ffs!
	FullMatrix<int, 6, 3> neighborBox, direction;
	FullMatrix<double, 6, 3> position;
	neighborBox << -1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1;
	direction = -neighborBox;
	position << 0., 0.5, 0.5, 1., 0.5, 0.5, 0.5, 0., 0.5, 0.5, 1., 0.5, 0.5, 0.5, 0., 0.5, 0.5, 1.;

	for (int iX = 0; iX < mDivisions[0]; ++iX)
		for (int iY = 0; iY < mDivisions[1]; ++iY)
			for (int iZ = 0; iZ < mDivisions[2]; ++iZ)
			{
				std::list<CollidableWallBase*> walls;
				SubBox& insideBox = *mSubBoxes[GetBoxIndex(iX, iY, iZ)];

				for (int iWall = 0; iWall < 6; ++iWall)
				{
					int indexX = neighborBox(iWall, 0);
					int indexY = neighborBox(iWall, 1);
					int indexZ = neighborBox(iWall, 2);

					SubBox& outsideBox = *mSubBoxes[GetBoxIndex(iX + indexX, iY + indexY, iZ + indexZ)];
					FullVector<double, 3> wallPosition( {
							(position(iWall, 0) + iX) * subBoxLength[0] + mSpecimen.GetBoundingBox()(0, 0),
							(position(iWall, 1) + iY) * subBoxLength[1] + mSpecimen.GetBoundingBox()(1, 0),
							(position(iWall, 2) + iZ) * subBoxLength[2] + mSpecimen.GetBoundingBox()(2, 0) });

					FullVector<double, 3> wallDirection;
					wallDirection << direction(iWall, 0), direction(iWall, 1), direction(iWall, 2);

					CollidableWallBase* wall;
					if (&insideBox == &outsideBox)
						wall = new CollidableWallPhysical(wallPosition, wallDirection, ++index);
					else
						wall = new CollidableWallVirtual(wallPosition, wallDirection, ++index);

					wall->SetBoxes(insideBox, outsideBox);
					walls.push_back(wall);

				}
				insideBox.SetWalls(walls);
			}
}

std::vector<NuTo::CollidableWallBase*> NuTo::SubBoxHandler::GetXYWalls(unsigned int rIndex)
{
	std::vector<CollidableWallBase*> xyWalls;

	for (auto wall : mSubBoxes[rIndex]->GetWalls())
		if (wall->GetDirection()[2] == 0)
			xyWalls.push_back(wall);

	if (xyWalls.size() != 4)
		throw NuTo::Exception("[NuTo::SubBoxHandler::BuildCylinder] This should not have happened. #BestThrow #JFU ");

	return xyWalls;
}

std::vector<NuTo::FullVector<double, Eigen::Dynamic> > NuTo::SubBoxHandler::GetXYCorners(
		std::vector<CollidableWallBase*> rWalls)
{
	double xMin = +INFINITY;
	double yMin = +INFINITY;
	double xMax = -INFINITY;
	double yMax = -INFINITY;
	double z = rWalls[0]->GetPosition()(2);

	for (auto wall : rWalls)
	{
		xMin = std::min(xMin, static_cast<double>(wall->GetPosition()(0)));
		yMin = std::min(yMin, static_cast<double>(wall->GetPosition()(1)));
		xMax = std::max(xMax, static_cast<double>(wall->GetPosition()(0)));
		yMax = std::max(yMax, static_cast<double>(wall->GetPosition()(1)));
	}

	std::vector<NuTo::FullVector<double, Eigen::Dynamic> > corners(4);

	corners[0] = FullVector<double, 3>( { xMin, yMin, z });
	corners[1] = FullVector<double, 3>( { xMin, yMax, z });
	corners[2] = FullVector<double, 3>( { xMax, yMin, z });
	corners[3] = FullVector<double, 3>( { xMax, yMax, z });

	return corners;
}

const std::vector<NuTo::SubBox*>& NuTo::SubBoxHandler::GetSubBoxes() const
{
	return mSubBoxes;
}

void NuTo::SubBoxHandler::AddSubBox(SubBox& rSubBox)
{
	mSubBoxes.push_back(&rSubBox);
}

std::vector<int> NuTo::SubBoxHandler::GetNInside(CollidableWallCylinder* rWallCylinder)
{
	std::vector<int> nInside(mSubBoxes.size());
	FullVector<double, 3> tmpNull( { 0., 0., 0. });
	for (unsigned int i = 0; i < mSubBoxes.size(); ++i)
	{
		std::vector<CollidableWallBase*> xyWalls = GetXYWalls(i);
		std::vector<NuTo::FullVector<double, Eigen::Dynamic> > xyCorners = GetXYCorners(xyWalls);
		nInside[i] = 0;
		for (int c = 0; c < 4; ++c)
		{
			CollidableParticleSphere sph(xyCorners[c], tmpNull, 0., 0., 0);
			if (rWallCylinder->IsInside(sph))
				nInside[i]++;
		}
	}
	return nInside;
}

void NuTo::SubBoxHandler::BuildCylinder(const int rNumThreads)
{
	mSubBoxes.clear();

	double radius = mSpecimen.GetLength(0) / 2.;
	double height = mSpecimen.GetLength(2);

	FullVector<double, 3> origin( {
			mSpecimen.GetBoundingBox()(0, 0) + radius,
			mSpecimen.GetBoundingBox()(1, 0) + radius,
			mSpecimen.GetBoundingBox()(2, 0) + height });

	FullVector<double, 3> direction( { 0., 0., 1 });

	CollidableWallCylinder* wallCylinder = new CollidableWallCylinder(origin, direction, radius, height, 1);

	// create box
	BuildSubBoxes(rNumThreads);

	AddSpheresToBoxes();

	// get number of corners inside the cylinder (nInside) for each cell
	std::vector<int> nInside = GetNInside(wallCylinder);

	// nInside = 1..3 --> add cylinder
	for (unsigned int i = 0; i < mSubBoxes.size(); ++i)
		if (1 <= nInside[i] && nInside[i] <= 3)
		{
			CollidableWallCylinder* subWallCylinder = new CollidableWallCylinder(origin, direction, radius, height, 1);
			subWallCylinder->SetBoxes(*mSubBoxes[i], *mSubBoxes[i]);
			mSubBoxes[i]->AddWall(*subWallCylinder);
		}

	// delete all physical side walls

	for (auto box : mSubBoxes)
	{
		std::vector<CollidableWallBase*> toDelete;
		for (auto wall : box->GetWalls())
		{
			bool isPhysicalWall = wall->IsPhysical();
			bool isSideWall = wall->GetDirection()(2) == 0;

			if (isPhysicalWall && isSideWall)
			{
				toDelete.push_back(wall);
			}
		}
		for (auto wall : toDelete)
		{
			box->RemoveWall(*wall);
			delete wall;
		}
	}
}

void NuTo::SubBoxHandler::AddSpheresToBoxes()
{

	double sTime = WallTime::Get();

	FullVector<double, 3> subLength = GetSubBoxLength();

	const int numParticles = mSpheres->GetNumParticles();

	std::cout << "[AddSpheresToBoxex]" << numParticles << std::endl;

	for (int sph = 0; sph < numParticles; ++sph)
	{
		std::vector<FullVector<double, 3> > corners(8, mSpheres->GetParticle(sph)->GetPosition());
		double radius = mSpheres->GetParticle(sph)->GetRadius();
		// add radius to get different points
		int cornerIndex = 0;
		for (int iX = -1; iX <= 1; iX += 2)
			for (int iY = -1; iY <= 1; iY += 2)
				for (int iZ = -1; iZ <= 1; iZ += 2)
				{
					corners[cornerIndex](0) += iX * radius;
					corners[cornerIndex](1) += iY * radius;
					corners[cornerIndex](2) += iZ * radius;
					cornerIndex++;
				}

		int xs = mDivisions[0];
		int ys = mDivisions[1];
		int zs = mDivisions[2];
		int xe = 0;
		int ye = 0;
		int ze = 0;

		for (int c = 0; c < 8; ++c)
		{
			// get xyz cell index of corner points
			int indX = std::floor(static_cast<double>((corners[c](0) - mSpecimen.GetBoundingBox()(0, 0)) / subLength[0]));
			int indY = std::floor(static_cast<double>((corners[c](1) - mSpecimen.GetBoundingBox()(1, 0)) / subLength[1]));
			int indZ = std::floor(static_cast<double>((corners[c](2) - mSpecimen.GetBoundingBox()(2, 0)) / subLength[2]));

			// get start and end indices of corner points
			xs = std::min(xs, indX);
			ys = std::min(ys, indY);
			zs = std::min(zs, indZ);
			xe = std::max(xe, indX);
			ye = std::max(ye, indY);
			ze = std::max(ze, indZ);

		}

		for (int indX = xs; indX <= xe; ++indX)
			for (int indY = ys; indY <= ye; ++indY)
				for (int indZ = zs; indZ <= ze; ++indZ)
				{

					// add sphere to corresponding box
					int boxIndex = mDivisions[1] * mDivisions[2] * indX + mDivisions[2] * indY + indZ;
					if (boxIndex < 0 || std::fabs(boxIndex) >= mSubBoxes.size())
						throw NuTo::Exception("[NuTo::SubBoxHandler::AddSpheresToBoxes] Box size/position does not match sphere size/position.");

					// only add new spheres
					CollidableBase* tmpColl = mSubBoxes[boxIndex]->GetCollidables().back();
					if (tmpColl != mSpheres->GetParticle(sph))
						mSubBoxes[boxIndex]->AddSphere(*mSpheres->GetParticle(sph));
				}

	}
	std::cout << "AddSpheresToBoxes() took " << WallTime::Get() - sTime
			<< " seconds." << std::endl;

}

void NuTo::SubBoxHandler::VisualizeBorders(std::string rFile)
{
#ifdef ENABLE_VISUALIZE
	NuTo::VisualizeUnstructuredGrid visuBorders;
	visuBorders.DefineCellDataVector("Direction");
	for (auto box : mSubBoxes)
		for (auto wall : box->GetWalls())
			wall->VisualizationStatic(visuBorders);

	visuBorders.ExportVtuDataFile(rFile);
#endif
}

double NuTo::SubBoxHandler::GetVolume() const
{
	return mSpecimen.GetVolume();
}

void NuTo::SubBoxHandler::PrintBoxes()
{
	for (auto box : mSubBoxes)
	{
		std::cout << "Box " << box->GetIndex() << "  =======================" << std::endl;
		box->Print();
	}

}
