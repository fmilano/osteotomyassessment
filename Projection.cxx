#include <limits>
#include <cmath>

#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkCellLocator.h>
#include <vtkMath.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkAppendPolyData.h>
#include <vtkSphereSource.h>
#include <vtkSTLWriter.h>

#include "vtkOBBTree2.h"

#include "Projection.h"

bool IsSxPieceSide(double piecePoint[], double planePoint[], double normal[], bool isNormalFlipped);

void Project(/*in*/ vtkPolyData* piece,
             /*in*/ vtkPolyData* plane,
             /*in*/ double normal[],
                    double planePush,
                    bool isNormalFlipped,
                    double threshold,
                    int projectionMode,
             /*out*/vtkPoints* intersectedPoints,
             /*out*/vtkIdList* intersectedCellsId,
             /*out*/vtkFloatArray* distances
            )
{

    vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();

    // Generates a tree to look for the cells
    vtkSmartPointer<vtkOBBTree2> tree = vtkSmartPointer<vtkOBBTree2>::New();
    tree->SetDataSet(piece);
    tree->BuildLocator();


    // Normalize the normal
    vtkMath::Normalize(normal);


    int notIntersected = 0;

    // generates a line segment that maybe intersects the piece
    double norm1[3];
    norm1[0] = normal[0] * (threshold + 1 + fabs(planePush));
    norm1[1] = normal[1] * (threshold + 1 + fabs(planePush));
    norm1[2] = normal[2] * (threshold + 1 + fabs(planePush));

    //intersectedPoints = vtkPoints::New();
    //vtkIdList* intersectedCells  = vtkIdList::New();

    // Iterates all the points
    // for each point in this grid, find the closest point in the piece
    for (vtkIdType i = 0; i < plane->GetNumberOfPoints(); i++) {
        double point0[3];
        plane->GetPoint(i, point0);

        double tmpPoint[3];
        tmpPoint[0] = point0[0] + normal[0] * planePush;
        tmpPoint[1] = point0[1] + normal[1] * planePush;
        tmpPoint[2] = point0[2] + normal[2] * planePush;

        if (i % 10 == 0)
        {
            vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
            sphere->SetCenter(tmpPoint);
            sphere->SetRadius(0.2);
            sphere->SetThetaResolution(8);
            sphere->SetPhiResolution(8);
            sphere->Update();

            appendPolyData->AddInput(sphere->GetOutput());
        }

        double a0[3];
        double a1[3];
        if (0 == projectionMode)
        {
          a0[0] = tmpPoint[0] + norm1[0];
          a0[1] = tmpPoint[1] + norm1[1];
          a0[2] = tmpPoint[2] + norm1[2];

          a1[0] = tmpPoint[0] - norm1[0];
          a1[1] = tmpPoint[1] - norm1[1];
          a1[2] = tmpPoint[2] - norm1[2];
        }
        else if (1 == projectionMode)
        {
          a0[0] = tmpPoint[0];
          a0[1] = tmpPoint[1];
          a0[2] = tmpPoint[2];

          a1[0] = tmpPoint[0] - norm1[0];
          a1[1] = tmpPoint[1] - norm1[1];
          a1[2] = tmpPoint[2] - norm1[2];
        }
        else if (2 == projectionMode)
        {
          a0[0] = tmpPoint[0] + norm1[0];
          a0[1] = tmpPoint[1] + norm1[1];
          a0[2] = tmpPoint[2] + norm1[2];

          a1[0] = tmpPoint[0];
          a1[1] = tmpPoint[1];
          a1[2] = tmpPoint[2];
        }

        vtkSmartPointer<vtkPoints> tmpIntersectedPoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkIdList> tmpIntersectedCells = vtkSmartPointer<vtkIdList>::New();

        tree->IntersectWithLine(a0, a1, tmpIntersectedPoints, tmpIntersectedCells);

        // detect the closest point in the piece
        double distance = std::numeric_limits<double>::max();
        double distance2 = std::numeric_limits<double>::max();
        double closestPoint[3];
        vtkIdType closestCellId;
        for (vtkIdType j = 0; j < tmpIntersectedPoints->GetNumberOfPoints(); j++) {
            double intersection[3];
            tmpIntersectedPoints->GetPoint(j, intersection);
            double tmpdist = std::sqrt(vtkMath::Distance2BetweenPoints(point0, intersection));
            double tmpdist2 = std::sqrt(vtkMath::Distance2BetweenPoints(tmpPoint, intersection));
            // Determine the side of the plane we are
            /*if (!IsSxPieceSide(intersection, point0, normal, isNormalFlipped)) {
                // We are not in the SxPiece side
                intersectedPoints->InsertNextPoint(intersection);
                intersectedCellsId->InsertNextId(tmpIntersectedCells->GetId(j));
                distances->InsertNextValue(tmpdist);
            }
            else*/
            if (tmpdist2 < distance2) {
                // save the closest point
                distance = tmpdist;
                distance2 = tmpdist2;
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
        if (IsSxPieceSide(closestPoint, point0, normal, isNormalFlipped))
            distances->InsertNextValue(-distance);
        else
            distances->InsertNextValue(distance);

        //colorimetryPoints->InsertNextPoint(point0);

        intersectedPoints->InsertNextPoint(closestPoint);
        intersectedCellsId->InsertNextId(closestCellId);
    }

    appendPolyData->Update();

    vtkSmartPointer<vtkSTLWriter> writer2 = vtkSmartPointer<vtkSTLWriter>::New();
    writer2->SetFileName("Prueba.stl");
    writer2->SetInput(appendPolyData->GetOutput());
    writer2->Write();

}

bool IsSxPieceSide(double piecePoint[], double planePoint[], double normal[], bool isNormalFlipped)
{
    // Sets a positive or negative distance according to the plane normal
    double cp[3];
    cp[0] = piecePoint[0] - planePoint[0];
    cp[1] = piecePoint[1] - planePoint[1];
    cp[2] = piecePoint[2] - planePoint[2];

    vtkMath::Normalize(cp);

    double tmpNormal[3];
    double factor = isNormalFlipped ? -1 : 1;
    tmpNormal[0] = normal[0] * factor;
    tmpNormal[1] = normal[1] * factor;
    tmpNormal[2] = normal[2] * factor;

   if (vtkMath::Dot(cp, tmpNormal) > 0)
       return true;

   return false;
}
