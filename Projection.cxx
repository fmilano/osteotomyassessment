#include <limits>
#include <cmath>

#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkCellLocator.h>
#include <vtkMath.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

#include "vtkOBBTree2.h"

#include "Projection.h"

bool IsSxPieceSide(double piecePoint[], double planePoint[], double normal[]);

void Project(/*in*/ vtkPolyData* piece,
             /*in*/ vtkPolyData* plane,
             /*in*/ double normal[],
                    double threshold,
             /*out*/vtkPoints* intersectedPoints,
             /*out*/vtkIdList* intersectedCellsId,
             /*out*/vtkFloatArray* distances
            )
{
    // Generates a tree to look for the cells
    vtkSmartPointer<vtkOBBTree2> tree = vtkSmartPointer<vtkOBBTree2>::New();
    tree->SetDataSet(piece);
    tree->BuildLocator();


    // Normalize the normal
    vtkMath::Normalize(normal);


    int notIntersected = 0;

    // generates a line segment that maybe intersects the piece
    double norm1[3];
    norm1[0] = normal[0] * (threshold + 1);
    norm1[1] = normal[1] * (threshold + 1);
    norm1[2] = normal[2] * (threshold + 1);

    //intersectedPoints = vtkPoints::New();
    //vtkIdList* intersectedCells  = vtkIdList::New();

    // Iterates all the points
    // for each point in this grid, find the closest point in the piece
    for (vtkIdType i = 0; i < plane->GetNumberOfPoints(); i++) {
        double point0[3];
        plane->GetPoint(i, point0);

        double a0[3] = {point0[0] + norm1[0],
                        point0[1] + norm1[1],
                        point0[2] + norm1[2]};
        double a1[3] = {point0[0] - norm1[0],
                        point0[1] - norm1[1],
                        point0[2] - norm1[2]};

        vtkSmartPointer<vtkPoints> tmpIntersectedPoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkIdList> tmpIntersectedCells = vtkSmartPointer<vtkIdList>::New();

        tree->IntersectWithLine(a0, a1, tmpIntersectedPoints, tmpIntersectedCells);

        // detect the closest point in the piece
        double distance = std::numeric_limits<double>::max();
        double closestPoint[3];
        vtkIdType closestCellId;
        for (vtkIdType j = 0; j < tmpIntersectedPoints->GetNumberOfPoints(); j++) {
            double intersection[3];
            tmpIntersectedPoints->GetPoint(j, intersection);
            double tmpdist = std::sqrt(vtkMath::Distance2BetweenPoints(point0, intersection));
            // Determine the side of the plane we are
            if (!IsSxPieceSide(intersection, point0, normal)) {
                // We are not in the SxPiece side
                intersectedPoints->InsertNextPoint(intersection);
                intersectedCellsId->InsertNextId(tmpIntersectedCells->GetId(j));
                distances->InsertNextValue(tmpdist);
            }
            else if (tmpdist < distance) {
                // save the closest point
                distance = tmpdist;
                closestPoint[0] = intersection[0];
                closestPoint[1] = intersection[1];
                closestPoint[2] = intersection[2];
                closestCellId = tmpIntersectedCells->GetId(j);
            }
        }

        // Points that did not intersect
        if (tmpIntersectedPoints->GetNumberOfPoints() == 0) {
            notIntersected++;
            continue;
        }

        // Discards the points outside the ROI
        if (distance > threshold)
            continue;

        // Sets a positive or negative distance according to the plane normal
        if (IsSxPieceSide(closestPoint, point0, normal))
            distances->InsertNextValue(-distance);
        else
            distances->InsertNextValue(distance);

        //colorimetryPoints->InsertNextPoint(point0);

        intersectedPoints->InsertNextPoint(closestPoint);
        intersectedCellsId->InsertNextId(closestCellId);
    }

}

bool IsSxPieceSide(double piecePoint[], double planePoint[], double normal[])
{
    // Sets a positive or negative distance according to the plane normal
    double cp[3];
    cp[0] = piecePoint[0] - planePoint[0];
    cp[1] = piecePoint[1] - planePoint[1];
    cp[2] = piecePoint[2] - planePoint[2];

    vtkMath::Normalize(cp);

   if (vtkMath::Dot(cp, normal) > 0)
       return true;

   return false;
}
