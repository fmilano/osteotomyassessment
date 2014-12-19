#ifndef PROJECTION_H
#define PROJECTION_H

#include <vtkPolyData.h>


void Project(/*in*/ vtkPolyData* piece,
             /*in*/ vtkPolyData* plane,
             /*in*/ double normal[],
                    double planePush,
                    bool isNormalFlipped,
                    double threshold,
                    int projectionMode,
            /*out*/ vtkPoints* intersectedPoints,
            /*out*/ vtkIdList* intersectedCellsId,
            /*out*/ vtkFloatArray* distances);

#endif
