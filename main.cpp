#include <iostream>
#include <cmath>
#include <limits>
#include <cstdlib>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkMetaImageWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPlaneSource.h>
#include <vtkCell.h>
#include <vtkTriangle.h>
#include <vtkPlane.h>
#include <vtkDoubleArray.h>
#include <vtkTriangleFilter.h>
#include <vtkCellLocator.h>
#include <vtkFloatArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyDataConnectivityFilter.h>


#include "vtkOBBTree2.h"
#include "vtkRANSACPlane.h"
#include "Projection.h"


vtkSmartPointer<vtkPlaneSource> GetPlaneParameters(vtkSmartPointer<vtkPolyData> inputPlane,
                                                   bool isInputPlane,
                                                   vtkSmartPointer<vtkCellLocator> tree,
                                                   double Spacing);
 
int main(int argc, char* argv[]) {

   if (argc < 4) {
      std::cerr << "Usage: " << std::endl << argv[0] << " specimen(.stl)"                    // argv[1]
                                                        " planned_osteotomy_virtual_plane(.stl)"   // argv[2]
                                                        " case_id"                           // argv[3]
                                                        " [inlier_threshold_in_mm)]"         // optional argv[4]
                                                        " [postSx_osteotomy_virtual_plane(.stl)]"  // optional argv[5]
                                                        " [results(.csv)]"                   // optional argv[6]
                                                        << std::endl;
      return 1;
   }

   std::string specimenFilename(argv[1]);
   std::string plannedOsteotomyPlaneFilename(argv[2]);
   std::string caseId(argv[3]);

   // The inlier threshold determines which points are considered to belong to the executed osteotomy plane
   double inlierThreshold = 0.5; // mm
   if (argc > 4)
   {
     inlierThreshold = atof(argv[4]);
   }

   // The postSxVirtualOsteotomy is a plane, usually set by hand, that is used to validate the algorithm
   std::string postSxOsteotomyPlaneFilename;
   if (argc > 5)
   {
     postSxOsteotomyPlaneFilename = (argv[5]);
   }

   // The results file contains the error distance for each sampled point
   std::string results("results_" + caseId + ".csv");
   if (argc > 6)
   {
      results = argv[6];
   }

   // RANSAC maximum interations
   int maxIterations = 5000;
   // Maximum distance considered when ray casting over the specimen (in mm)
   double regionOfInterestThreshold = 20;
   // Spacing of the plane resampling (in mm)
   double spacing = 0.1;


   std::cout << "**************************************************************"
                "******************" << std::endl;
   std::cout << "Specimen: " << specimenFilename << " Planned Osteotomy Virtual Plane: " << plannedOsteotomyVirtualPlaneFilename << std::endl;
   std::cout << "Inlier threshold: " << inlierThreshold << std::endl;

   std::cout << "--------------------------------------------------------------"
                "------------------" << std::endl;

   vtkSmartPointer<vtkSTLReader> specimenReader = vtkSmartPointer<vtkSTLReader>::New();
   specimenReader->SetFileName(specimenFilename);
   specimenReader->Update();

   vtkSmartPointer<vtkPolyData> specimen = pieceReader->GetOutput();

   vtkSmartPointer<vtkSTLReader> plannedOsteotomyPlaneReader = vtkSmartPointer<vtkSTLReader>::New();
   plannedOsteotomyPlaneReader->SetFileName(plannedOsteotomyPlaneFilename);
   plannedOsteotomyPlaneReader->Update();

  vtkSmartPointer<vtkPolyData> plannedOsteotomyPlane = plannedOsteotomyPlaneReader->GetOutput();

  /////////////////////////////////////////////////////////////////////////////
  // Calculates plane parameters
  /////////////////////////////////////////////////////////////////////////////

   // Generates a tree to look for the cells
   vtkSmartPointer<vtkOBBTree2> tree = vtkSmartPointer<vtkOBBTree2>::New();
   tree->SetDataSet(specimen);
   tree->BuildLocator();

   vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
   cellLocator->SetDataSet(specimen);
   cellLocator->BuildLocator();

   // The planned osteotomy plane is a thick plane, so we need a cellLocator to localize the cell
   // of the specimen that is closest to the plane; then we can return the closest face of the plane.
   vtkSmartPointer<vtkPlaneSource> planeSource = GetPlaneParameters(plannedOsteotomyPlane,
                                                                    true,
                                                                    cellLocator,
                                                                    spacing);

   vtkSmartPointer<vtkPolyData> planeData = planeSource->GetOutput();

   /////////////////////////////////////////////////////////////////////////////
   // Project points
   /////////////////////////////////////////////////////////////////////////////
   vtkSmartPointer<vtkPoints> intersectedPoints = vtkSmartPointer<vtkPoints>::New();
   vtkSmartPointer<vtkIdList> intersectedCellsId = vtkSmartPointer<vtkIdList>::New();

   vtkSmartPointer<vtkFloatArray> distances = vtkSmartPointer<vtkFloatArray>::New();
   distances->SetNumberOfComponents(1);
   distances->SetName("Distances");

   Project(specimen, planeData, planeSource->GetNormal(), regionOfInterestThreshold, intersectedPoints, intersectedCellsId, distances);

   /////////////////////////////////////////////////////////////////////////////
   // Prune with RANSAC
   /////////////////////////////////////////////////////////////////////////////
   vtkSmartPointer<vtkPolyData> executedOsteotomyPointCloud = vtkSmartPointer<vtkPolyData>::New();
   executedOsteotomyPointCloud->SetPoints(intersectedPoints);
   executedOsteotomyPointCloud->Update();

   vtkSmartPointer<vtkRANSACPlane> ransacPlane = vtkSmartPointer<vtkRANSACPlane>::New();
   ransacPlane->SetInput(executedOsteotomyPointCloud);

   double bounds[6];
   planeData->ComputeBounds();
   planeData->GetBounds(bounds);

   ransacPlane->SetBounds(bounds);
   ransacPlane->SetInlierThreshold(inlierThreshold);
   ransacPlane->SetMaxIterations(maxIterations);
   ransacPlane->Update();

   /////////////////////////////////////////////////////////////////////////////
   // Create colorimetry
   /////////////////////////////////////////////////////////////////////////////
   vtkSmartPointer<vtkPolyData> colorimetry = vtkSmartPointer<vtkPolyData>::New();
   colorimetry->ShallowCopy(inputPiece);

   // Create color data
   vtkSmartPointer<vtkFloatArray> colorData = vtkSmartPointer<vtkFloatArray>::New();
   colorData->SetNumberOfComponents(1);
   colorData->SetName("Colors");
   colorData->SetNumberOfTuples(inputPiece->GetNumberOfCells());

   for (vtkIdType i = 0; i < inputPiece->GetNumberOfCells(); i++)
       colorData->SetValue(i, std::numeric_limits<float>::quiet_NaN());


   // Paint all the positive distances (non specimen side); is helps visualization
   for (vtkIdType i = 0; i < distances->GetNumberOfTuples(); i++) {
        float distance = distances->GetValue(i);
        if (distance > 0)
            colorData->SetValue(intersectedCellsId->GetId(i), distance);
   }

   std::vector<unsigned int> maxInlierIndices = ransacPlane->MaxInlierIndices;

   // Color just the RANSAC computed indices
   for (std::vector<unsigned int>::iterator it = maxInlierIndices.begin(); it != maxInlierIndices.end(); it++) {
       colorData->SetValue(intersectedCellsId->GetId(*it), distances->GetValue(*it));
   }

   colorimetry->GetCellData()->SetScalars(colorData);
      
   /////////////////////////////////////////////////////////////////////////////
   // Results
   /////////////////////////////////////////////////////////////////////////////
   std::cout << "Results:" << std::endl;

   vtkSmartPointer<vtkPlane> bestPlane = ransacPlane->GetBestPlane();

   double nInput[3], nRansac[3], nManual[3];
   bestPlane->GetNormal(nRansac);
   //manualPlane->GetNormal(nManual);
   planeSource->GetNormal(nInput);

   vtkMath::Normalize(nRansac);
   //vtkMath::Normalize(nManual);
   vtkMath::Normalize(nInput);

   /*vtkSmartPointer<vtkSTLReader> manualReader = vtkSmartPointer<vtkSTLReader>::New();
   manualReader->SetFileName(argv[5]);
   manualReader->Update();

  vtkSmartPointer<vtkPolyData> inputManualPlane = manualReader->GetOutput();
 */

  /* double angleInputRansac = acos(vtkMath::Dot(nInput, nRansac))*180.0/M_PI;
   angleInputRansac = (angleInputRansac >= 90 ? 180 - angleInputRansac : angleInputRansac);
   double angleInputManual = acos(vtkMath::Dot(nInput, nManual))*180.0/M_PI;
   angleInputManual = (angleInputManual >= 90 ? 180 - angleInputManual : angleInputManual);
   std::cout << "Angle input - ransac: " << angleInputRansac << " degrees." << std::endl;
   std::cout << "Angle input - manual: " << angleInputManual << " degrees." << std::endl;
*/
   // calculate distance
   double cInput[3], cRansac[3], cManual[3]; // projInputRansac[3], projInputManual[3];
   bestPlane->GetOrigin(cRansac);
   //manualPlane->GetCenter(cManual);
   planeSource->GetCenter(cInput);



  /* double distanceRansacInput = vtkPlane::DistanceToPlane(cInput, nRansac, cRansac);
   double distanceManualInput = vtkPlane::DistanceToPlane(cInput, nManual, cManual);

   std::cout << "Distance input - ransac: " << distanceRansacInput << " mm" << std::endl;
   std::cout << "Distance input - manual: " << distanceManualInput << " mm" << std::endl;

   fstream results;

   results.open(argv[6], fstream::out | fstream::app);


   results << argv[3] << ", "
           << InlierThreshold  << ","
           << angleInputManual << ","
           << angleInputRansac << ","
           << distanceManualInput << ","
           << distanceRansacInput << std::endl;

   results.close();
*/
   std::cout << "**************************************************************"
                "******************" << std::endl;

   vtkSmartPointer<vtkTriangleFilter> tf = vtkSmartPointer<vtkTriangleFilter>::New();
   tf->SetInput(ransacPlane->GetOutput());
   tf->Update();

   vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
//writer->SetFileName(std::string(argv[3]).append("_found.stl").c_str());
   std::string foundFileName(argv[2]);
   // remove extension
   writer->SetFileName(foundFileName.erase(foundFileName.length() - 4).append("_found.stl").c_str());
   writer->SetInput(tf->GetOutput());
   writer->Write();

   tf->SetInput(planeData);
   tf->Update();

   std::string singleSideFileName(argv[2]);
   vtkSmartPointer<vtkSTLWriter> writer2 = vtkSmartPointer<vtkSTLWriter>::New();
   writer2->SetFileName(singleSideFileName.erase(singleSideFileName.length() - 4).append("_single_side.stl").c_str());
   writer2->SetInput(tf->GetOutput());
   writer2->Write();

   std::string colorimetryFileName(argv[2]);
   vtkSmartPointer<vtkXMLPolyDataWriter> writer3 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
   writer3->SetFileName(colorimetryFileName.erase(colorimetryFileName.length() - 4).append("_colorimetry.vtp").c_str());
   writer3->SetInput(colorimetry);
   writer3->Write();

   return 0;
}

/////////////////////////////////////////////////////////////////////////////
// Extract one side of the plane
/////////////////////////////////////////////////////////////////////////////
vtkSmartPointer<vtkPlaneSource> GetPlaneParameters(vtkSmartPointer<vtkPolyData> inputPlane,
                                                   bool isThick,
                                                   vtkSmartPointer<vtkCellLocator> tree,
                                                   double spacing) {

   double maxArea = 0;
   double maxp0[3], maxp1[3], maxp2[3];
   for (vtkIdType i = 0; i < inputPlane->GetNumberOfCells(); i++) {
      vtkCell* cell = inputPlane->GetCell(i);

      vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(cell);
      double p0[3], p1[3], p2[3];
      triangle->GetPoints()->GetPoint(0, p0);
      triangle->GetPoints()->GetPoint(1, p1);
      triangle->GetPoints()->GetPoint(2, p2);

      double t = vtkTriangle::TriangleArea(p0, p1, p2);
      if (t > maxArea) {
         maxArea = t;

         maxp0[0] = p0[0];
         maxp0[1] = p0[1];
         maxp0[2] = p0[2];

         maxp1[0] = p1[0];
         maxp1[1] = p1[1];
         maxp1[2] = p1[2];

         maxp2[0] = p2[0];
         maxp2[1] = p2[1];
         maxp2[2] = p2[2];
      }
   }

   double distp0p1 = std::sqrt(vtkMath::Distance2BetweenPoints(maxp0, maxp1));
   double distp0p2 = std::sqrt(vtkMath::Distance2BetweenPoints(maxp0, maxp2));
   double distp1p2 = std::sqrt(vtkMath::Distance2BetweenPoints(maxp1, maxp2));

   double origin[3], lateral1[3], lateral2[3];

   double maxDist = distp1p2;

   origin[0]   = maxp0[0];
   origin[1]   = maxp0[1];
   origin[2]   = maxp0[2];

   lateral1[0] = maxp1[0];
   lateral1[1] = maxp1[1];
   lateral1[2] = maxp1[2];

   lateral2[0] = maxp2[0];
   lateral2[1] = maxp2[1];
   lateral2[2] = maxp2[2];

   if (distp0p1 > maxDist) {
      origin[0]   = maxp2[0];
      origin[1]   = maxp2[1];
      origin[2]   = maxp2[2];

      lateral2[0] = maxp0[0];
      lateral2[1] = maxp0[1];
      lateral2[2] = maxp0[2];

      maxDist = distp0p1;
   }
   if (distp0p2 > maxDist) {
      origin[0]   = maxp1[0];
      origin[1]   = maxp1[1];
      origin[2]   = maxp1[2];

      lateral1[0] = maxp0[0];
      lateral1[1] = maxp0[1];
      lateral1[2] = maxp0[2];

      lateral2[0] = maxp2[0];
      lateral2[1] = maxp2[1];
      lateral2[2] = maxp2[2];

      maxDist = distp0p2;
   }

   // create the plane

   double normal[3];
   vtkTriangle::ComputeNormalDirection(origin, lateral1, lateral2, normal);

   vtkMath::Normalize(normal);

   vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();

   double projOrigin[3], projLateral1[3], projLateral2[3];

   double center[3];
   inputPlane->GetCenter(center);

   vtkPlane::ProjectPoint(origin,   center, normal, projOrigin);
   vtkPlane::ProjectPoint(lateral1, center, normal, projLateral1);
   vtkPlane::ProjectPoint(lateral2, center, normal, projLateral2);

   planeSource->SetOrigin(projOrigin);
   planeSource->SetPoint1(projLateral1);
   planeSource->SetPoint2(projLateral2);

   planeSource->GetNormal(normal);

   if (isInputPlane) {
      // At this point, the plane is centered in the middle of the box
      // We create two points to detect automatically if we have to flip
      // the sides
      double pointFlip[3], pointNotFlip[3], closestPoint[3];
      pointFlip[0] = center[0] + normal[0] * -20;
      pointFlip[1] = center[1] + normal[1] * -20;
      pointFlip[2] = center[2] + normal[2] * -20;

      pointNotFlip[0] = center[0] + normal[0] * 20;
      pointNotFlip[1] = center[1] + normal[1] * 20;
      pointNotFlip[2] = center[2] + normal[2] * 20;

      double distFlip, distNotFlip;
      int subId;
      vtkIdType idType;

      tree->FindClosestPoint(pointFlip, closestPoint, idType, subId, distFlip);
      tree->FindClosestPoint(pointNotFlip, closestPoint, idType, subId, distNotFlip);
      std::cout << "distFlip: " << distFlip << "; distNoFlip:" << distNotFlip << std::endl;
      if (distFlip < distNotFlip) {
         normal[0] = -normal[0];
         normal[1] = -normal[1];
         normal[2] = -normal[2];

         planeSource->SetNormal(normal);
         planeSource->Update();
      }

      planeSource->Push(1);

     // calculates the resolution
     int xresolution = std::sqrt(vtkMath::Distance2BetweenPoints(projOrigin, projLateral1)) / Spacing;
     int yresolution = std::sqrt(vtkMath::Distance2BetweenPoints(projOrigin, projLateral2)) / Spacing;

     planeSource->SetResolution(xresolution, yresolution);
   }
   planeSource->Update();

   return planeSource;
}
