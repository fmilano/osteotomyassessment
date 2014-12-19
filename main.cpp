#include <iostream>
#include <algorithm>
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
#include <vtkAppendPolyData.h>
#include <vtkSphereSource.h>


#include "vtkOBBTree2.h"
#include "vtkLORANSACPlane.h"
#include "Projection.h"

struct PlaneParameters
{
  vtkSmartPointer<vtkPlaneSource> planeSource;
  bool isNormalFlipped;
};


PlaneParameters GetPlaneParameters(vtkSmartPointer<vtkPolyData> inputPlane,
                                                   bool isInputPlane = false,
                                                   bool forceSideFlip = false,
                                                   vtkSmartPointer<vtkCellLocator> tree = 0,
                                                   double Spacing = 0);
 
int main(int argc, char* argv[]) {

   if (argc < 4) {
      std::cerr << "Usage: " << std::endl << argv[0] << " specimen(.stl)"                         // argv[1]
                                                        " planned_osteotomy_virtual_plane(.stl)"  // argv[2]
                                                        " case_id"                                // argv[3]
                                                        " projection_mode"                        // argv[4]
                                                        " [forceSideFlip]"                        // argv[5]
                                                        " [spacing_in_mm]"                        // optional argv[6]
                                                        " [plane_push_in_mm]"                     // optional argv[7]
                                                        " [ROI_threshold]"                        // optional argv[8]
                                                        " [inlier_threshold_in_mm]"               // optional argv[9]
                                                        " [postSx_osteotomy_virtual_plane(.stl)]" // optional argv[10]
                                                        << std::endl;
      return 1;
   }

   std::string specimenFilename(argv[1]);
   std::string plannedOsteotomyPlaneFilename(argv[2]);
   std::string caseId(argv[3]);
   int projectionMode = atoi(argv[4]);

   bool forceSideFlip = false;
   if (argc > 5)
   {
     if (atoi(argv[5]) > 0)
       forceSideFlip = true;
   }

   // Spacing of the plane resampling
   double spacing = 0.1; // mm
   if (argc > 6)
   {
     spacing = atof(argv[6]);
   }   

   // Push for plane projection
   double planePush = 0; // mm
   if (argc > 7)
   {
     planePush = atof(argv[7]);
   }

   // Maximum distance considered when ray casting over the specimen (in mm)
   double regionOfInterestThreshold = 20;
   if (argc > 8)
   {
     regionOfInterestThreshold = atof(argv[8]);
   }

   // The inlier threshold determines which points are considered to belong to the executed osteotomy plane
   double inlierThreshold = 0.5; // mm
   if (argc > 9)
   {
     inlierThreshold = atof(argv[9]);
   }

   // The postSxVirtualOsteotomy is a plane, usually set by hand, that is used to validate the algorithm
   std::string postSxOsteotomyPlaneFilename;
   if (argc > 10)
   {
     postSxOsteotomyPlaneFilename = (argv[10]);
   }

   // RANSAC maximum interations
   int maxIterations = 5000;

   std::cout << "--------------------------------------------------------------"
                "------------------" << std::endl;
   std::cout << "Case: " << caseId << std::endl;
   std::cout << "--------------------------------------------------------------"
                "------------------" << std::endl;
   std::cout << "Specimen: " << specimenFilename << std::endl;
   std::cout << "Planned Osteotomy Virtual Plane: " << plannedOsteotomyPlaneFilename << std::endl;
   std::cout << "Spacing: " << spacing << "mm" << std::endl;
   std::cout << "Inlier threshold: " << inlierThreshold << "mm" << std::endl;

   vtkSmartPointer<vtkSTLReader> specimenReader = vtkSmartPointer<vtkSTLReader>::New();
   specimenReader->SetFileName(specimenFilename.c_str());
   specimenReader->Update();

   vtkSmartPointer<vtkPolyData> specimen = specimenReader->GetOutput();

   vtkSmartPointer<vtkSTLReader> plannedOsteotomyPlaneReader = vtkSmartPointer<vtkSTLReader>::New();
   plannedOsteotomyPlaneReader->SetFileName(plannedOsteotomyPlaneFilename.c_str());
   plannedOsteotomyPlaneReader->Update();

  vtkSmartPointer<vtkPolyData> plannedOsteotomyPlane = plannedOsteotomyPlaneReader->GetOutput();

  /////////////////////////////////////////////////////////////////////////////
  // Calculates plane parameters
  /////////////////////////////////////////////////////////////////////////////
  std::cout << "Processing input data";

   // Generates a tree to look for the cells
   vtkSmartPointer<vtkOBBTree2> tree = vtkSmartPointer<vtkOBBTree2>::New();
   tree->SetDataSet(specimen);
   tree->BuildLocator();

   vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
   cellLocator->SetDataSet(specimen);
   cellLocator->BuildLocator();

   // The planned osteotomy plane is a thick plane, so we need a cellLocator to localize the cell
   // of the specimen that is closest to the plane; then we can return the closest face of the plane.
   PlaneParameters planeParameters = GetPlaneParameters(plannedOsteotomyPlane,
                                                        true,
                                                        forceSideFlip,
                                                        cellLocator,
                                                        spacing);

   vtkSmartPointer<vtkPlaneSource> planeSource = planeParameters.planeSource;
   bool isNormalFlipped = planeParameters.isNormalFlipped;

   vtkSmartPointer<vtkPolyData> planeData = planeSource->GetOutput();

   std::cout << " [ok]" << std::endl;
   /////////////////////////////////////////////////////////////////////////////
   // Project points
   /////////////////////////////////////////////////////////////////////////////
   std::cout << "Projecting points";

   vtkSmartPointer<vtkPoints> intersectedPoints = vtkSmartPointer<vtkPoints>::New();
   vtkSmartPointer<vtkIdList> intersectedCellsId = vtkSmartPointer<vtkIdList>::New();

   vtkSmartPointer<vtkFloatArray> distances = vtkSmartPointer<vtkFloatArray>::New();
   distances->SetNumberOfComponents(1);
   distances->SetName("Distances");

   Project(specimen, planeData, planeSource->GetNormal(), planePush, isNormalFlipped, regionOfInterestThreshold, projectionMode, intersectedPoints, intersectedCellsId, distances);

   std::cout << " [ok]" << std::endl;
   /////////////////////////////////////////////////////////////////////////////
   // Prune with RANSAC
   /////////////////////////////////////////////////////////////////////////////
   std::cout << "Prunning data with LO-RANSAC";

   vtkSmartPointer<vtkPolyData> executedOsteotomyPointCloud = vtkSmartPointer<vtkPolyData>::New();
   executedOsteotomyPointCloud->SetPoints(intersectedPoints);
   executedOsteotomyPointCloud->Update();

   vtkSmartPointer<vtkLORANSACPlane> ransacPlane = vtkSmartPointer<vtkLORANSACPlane>::New();
   ransacPlane->SetInput(executedOsteotomyPointCloud);

   double bounds[6];
   planeData->ComputeBounds();
   planeData->GetBounds(bounds);

   ransacPlane->SetBounds(bounds);
   ransacPlane->SetInlierThreshold(inlierThreshold);
   ransacPlane->SetMaxIterations(maxIterations);
   ransacPlane->SetGoodEnough(0.95);
   ransacPlane->Update();

   std::vector<unsigned int> maxInlierIndices = ransacPlane->MaxInlierIndices;

   std::cout << " [ok]" << std::endl;
   std::cout << "Results:" << std::endl;


   std::string planeName = plannedOsteotomyPlaneFilename.erase(plannedOsteotomyPlaneFilename.length() - 4);
   /////////////////////////////////////////////////////////////////////////////
   // Create SPECIMEN big points STL
   /////////////////////////////////////////////////////////////////////////////
   std::cout << "    Saving SPECIMEN big points STL";

   std::vector<unsigned int> displayedPointsIndices;
   vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();

   for (std::vector<unsigned int>::iterator it = maxInlierIndices.begin(); it != maxInlierIndices.end(); it++)
   {
     double center[3];
     intersectedPoints->GetPoint(*it, center);
     bool shouldDisplayPoint = true;
     for (std::vector<unsigned int>::iterator it2 = displayedPointsIndices.begin(); it2 != displayedPointsIndices.end(); it2++)
     {
       double center2[3];
       intersectedPoints->GetPoint(*it2, center2);

       if (sqrt(vtkMath::Distance2BetweenPoints(center, center2)) < 3)
         shouldDisplayPoint = false;
     }

     if (shouldDisplayPoint)
     {
       vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
       sphere->SetCenter(center);
       sphere->SetRadius(0.8);
       sphere->SetThetaResolution(16);
       sphere->SetPhiResolution(16);
       sphere->Update();

       appendPolyData->AddInput(sphere->GetOutput());
       displayedPointsIndices.push_back(*it);
     }
   }

   appendPolyData->Update();

   vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
   writer->SetFileName((caseId + "-" + planeName + "_SpecimenBigPts.stl").c_str());
   writer->SetInput(appendPolyData->GetOutput());
   writer->Write();

   std::cout << " [ok]" << std::endl;

   /////////////////////////////////////////////////////////////////////////////
   // Create PLANE big points STL
   /////////////////////////////////////////////////////////////////////////////
   std::cout << "    Saving PLANE big points STL";

   std::vector<unsigned int> displayedPointsIndices2;
   vtkSmartPointer<vtkAppendPolyData> appendPolyData2 = vtkSmartPointer<vtkAppendPolyData>::New();

   for (vtkIdType i = 0; i < planeData->GetNumberOfPoints(); ++i)
   {
     double center[3];
     planeData->GetPoint(i, center);
     bool shouldDisplayPoint = true;
     for (std::vector<unsigned int>::iterator it2 = displayedPointsIndices2.begin(); it2 != displayedPointsIndices2.end(); it2++)
     {
       double center2[3];
       planeData->GetPoint(*it2, center2);

       if (sqrt(vtkMath::Distance2BetweenPoints(center, center2)) < 8)
         shouldDisplayPoint = false;
     }

     if (shouldDisplayPoint)
     {
       vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
       sphere->SetCenter(center);
       sphere->SetRadius(0.8);
       sphere->SetThetaResolution(16);
       sphere->SetPhiResolution(16);
       sphere->Update();

       appendPolyData2->AddInput(sphere->GetOutput());
       displayedPointsIndices2.push_back(i);
     }
   }

   appendPolyData2->Update();

   vtkSmartPointer<vtkSTLWriter> writer2 = vtkSmartPointer<vtkSTLWriter>::New();
   writer2->SetFileName((caseId + "-" + planeName + "_PlaneBigPts.stl").c_str());
   writer2->SetInput(appendPolyData2->GetOutput());
   writer2->Write();

   std::cout << " [ok]" << std::endl;

   /////////////////////////////////////////////////////////////////////////////
   // Create colorimetry
   /////////////////////////////////////////////////////////////////////////////
   std::cout << "    Saving colorimetry";

   vtkSmartPointer<vtkPolyData> colorimetry = vtkSmartPointer<vtkPolyData>::New();
   colorimetry->ShallowCopy(specimen);

   // Create color data
   vtkSmartPointer<vtkFloatArray> colorData = vtkSmartPointer<vtkFloatArray>::New();
   colorData->SetNumberOfComponents(1);
   colorData->SetName("Colors");
   colorData->SetNumberOfTuples(specimen->GetNumberOfCells());

   for (vtkIdType i = 0; i < specimen->GetNumberOfCells(); i++)
       colorData->SetValue(i, std::numeric_limits<float>::quiet_NaN());

   // Paint all the positive distances, even thou they are not pruned by RANSAC; it helps visualization
   /*for (vtkIdType i = 0; i < distances->GetNumberOfTuples(); i++) {
        float distance = distances->GetValue(i);
        if (distance > 0)
            colorData->SetValue(intersectedCellsId->GetId(i), distance);

   }
   */

   // Color just the RANSAC computed indices
   for (std::vector<unsigned int>::iterator it = maxInlierIndices.begin(); it != maxInlierIndices.end(); it++) {
       colorData->SetValue(intersectedCellsId->GetId(*it), distances->GetValue(*it));
   }

   colorimetry->GetCellData()->SetScalars(colorData);

   vtkSmartPointer<vtkXMLPolyDataWriter> xmlWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
   xmlWriter->SetFileName((caseId + "-" + planeName + "_colorimetry.vtp").c_str());
   xmlWriter->SetInput(colorimetry);
   xmlWriter->Write();

   std::cout << " [ok]" << std::endl;
      
   /////////////////////////////////////////////////////////////////////////////
   // Save result distances in CSV file
   /////////////////////////////////////////////////////////////////////////////

   std::cout << "    Saving distances";

   ofstream distancesResults;
   distancesResults.open((caseId + "-" + planeName + "_distances.csv").c_str(), std::ofstream::trunc);
   for (std::vector<unsigned int>::iterator it = maxInlierIndices.begin(); it != maxInlierIndices.end(); it++)
   {
     distancesResults << distances->GetValue(*it) << std::endl;
   }
   distancesResults.close();

   std::cout << " [ok]" << std::endl;

   /////////////////////////////////////////////////////////////////////////////
   // Save command line
   /////////////////////////////////////////////////////////////////////////////

   std::cout << "    Saving command line";

   ofstream commandLine;
   commandLine.open((caseId + "-" + planeName + ".cmd").c_str(), std::ofstream::trunc);
   for (int i = 0; i < argc; i++)
   {
       commandLine << argv[i] << " ";
   }
   commandLine.close();

   std::cout << " [ok]" << std::endl;

   /////////////////////////////////////////////////////////////////////////////
   // Save validation data
   /////////////////////////////////////////////////////////////////////////////

   if (!postSxOsteotomyPlaneFilename.empty())
   {
     std::cout << "    Saving validation results";
     vtkSmartPointer<vtkPlane> bestPlane = ransacPlane->GetBestPlane();

     vtkSmartPointer<vtkSTLReader> postSxOsteotomyPlaneReader = vtkSmartPointer<vtkSTLReader>::New();
     postSxOsteotomyPlaneReader->SetFileName(postSxOsteotomyPlaneFilename.c_str());
     postSxOsteotomyPlaneReader->Update();

     vtkSmartPointer<vtkPolyData> postSxOsteotomyPlane = postSxOsteotomyPlaneReader->GetOutput();
     vtkSmartPointer<vtkPlaneSource> manualPlane = GetPlaneParameters(postSxOsteotomyPlane).planeSource;

     double nInput[3], nRansac[3], nManual[3];
     bestPlane->GetNormal(nRansac);
     manualPlane->GetNormal(nManual);
     planeSource->GetNormal(nInput);

     vtkMath::Normalize(nRansac);
     vtkMath::Normalize(nManual);
     vtkMath::Normalize(nInput);

     double angleInputRansac = acos(vtkMath::Dot(nInput, nRansac))*180.0/M_PI;
     angleInputRansac = (angleInputRansac >= 90 ? 180 - angleInputRansac : angleInputRansac);
     double angleInputManual = acos(vtkMath::Dot(nInput, nManual))*180.0/M_PI;
     angleInputManual = (angleInputManual >= 90 ? 180 - angleInputManual : angleInputManual);
     //std::cout << "Angle input - ransac: " << angleInputRansac << " degrees." << std::endl;
     //std::cout << "Angle input - manual: " << angleInputManual << " degrees." << std::endl;

     // calculate distance
     double cInput[3], cRansac[3], cManual[3]; // projInputRansac[3], projInputManual[3];
     bestPlane->GetOrigin(cRansac);
     manualPlane->GetCenter(cManual);
     planeSource->GetCenter(cInput);

     double distanceRansacInput = vtkPlane::DistanceToPlane(cInput, nRansac, cRansac);
     double distanceManualInput = vtkPlane::DistanceToPlane(cInput, nManual, cManual);

     std::cout << "Distance input - ransac: " << distanceRansacInput << " mm" << std::endl;
     std::cout << "Distance input - manual: " << distanceManualInput << " mm" << std::endl;

     ofstream validationResults;
     validationResults.open((caseId + "-" + planeName +  "_validation.csv").c_str(), std::ofstream::trunc);
     validationResults << caseId   << ", "
             << spacing  << ","
             << inlierThreshold  << ","
             << angleInputManual << ","
             << angleInputRansac << ","
             << distanceManualInput << ","
             << distanceRansacInput << std::endl;

     validationResults.close();
     std::cout << " [ok]" << std::endl;
   }

   /////////////////////////////////////////////////////////////////////////////
   // Save found plane
   /////////////////////////////////////////////////////////////////////////////
   std::cout << "    Saving found plane";

   double foundPlaneNormal[3], foundPlaneCenter[3], projOrigin[3], projLateral1[3], projLateral2[3];
   ransacPlane->GetBestPlane()->GetNormal(foundPlaneNormal);
   ransacPlane->GetBestPlane()->GetOrigin(foundPlaneCenter);

   vtkPlane::ProjectPoint(planeSource->GetOrigin(), foundPlaneCenter, foundPlaneNormal, projOrigin);
   vtkPlane::ProjectPoint(planeSource->GetPoint1(), foundPlaneCenter, foundPlaneNormal, projLateral1);
   vtkPlane::ProjectPoint(planeSource->GetPoint2(), foundPlaneCenter, foundPlaneNormal, projLateral2);

   vtkSmartPointer<vtkPlaneSource> foundPlane = vtkSmartPointer<vtkPlaneSource>::New();
   foundPlane->SetOrigin(projOrigin);
   foundPlane->SetPoint1(projLateral1);
   foundPlane->SetPoint2(projLateral2);
   foundPlane->Update();

   vtkSmartPointer<vtkTriangleFilter> tf = vtkSmartPointer<vtkTriangleFilter>::New();
   tf->SetInput(foundPlane->GetOutput());
   tf->Update();

   writer = vtkSmartPointer<vtkSTLWriter>::New();
   writer->SetFileName((caseId + "-" + planeName + "_found.stl").c_str());
   writer->SetInput(tf->GetOutput());
   writer->Write();

   std::cout << " [ok]" << std::endl;

   /////////////////////////////////////////////////////////////////////////////
   // Save single side plane used for computing the cut plane
   /////////////////////////////////////////////////////////////////////////////
   std::cout << "    Saving single side of planned plane";

   tf = vtkSmartPointer<vtkTriangleFilter>::New();
   tf->SetInput(planeData);
   tf->Update();

   writer = vtkSmartPointer<vtkSTLWriter>::New();
   writer->SetFileName((caseId + "-" + planeName + "_single_side.stl").c_str());
   writer->SetInput(tf->GetOutput());
   writer->Write();

   std::cout << " [ok]" << std::endl;
   std::cout << "               <<<<<<<>>>>>>>" << std::endl;

   return 0;
}

/////////////////////////////////////////////////////////////////////////////
// Extract one side of the plane
/////////////////////////////////////////////////////////////////////////////
PlaneParameters GetPlaneParameters(vtkSmartPointer<vtkPolyData> inputPlane,
                                                   bool isThick,
                                                   bool forceSideFlip,
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

   distp0p1 = std::sqrt(vtkMath::Distance2BetweenPoints(origin, lateral1));
   distp0p2 = std::sqrt(vtkMath::Distance2BetweenPoints(origin, lateral2));

   if (distp0p2 > distp0p1)
   {
     std::swap(lateral1[0], lateral2[0]);
     std::swap(lateral1[1], lateral2[1]);
     std::swap(lateral1[2], lateral2[2]);
   }

   // create the plane

   vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
   planeSource->SetOrigin(origin);
   planeSource->SetPoint1(lateral1);
   planeSource->SetPoint2(lateral2);

   double normal[3];
   planeSource->GetNormal(normal);
   vtkMath::Normalize(normal);

   double projOrigin[3], projLateral1[3], projLateral2[3];
   double center[3];
   inputPlane->GetCenter(center);

   vtkPlane::ProjectPoint(origin,   center, normal, projOrigin);
   vtkPlane::ProjectPoint(lateral1, center, normal, projLateral1);
   vtkPlane::ProjectPoint(lateral2, center, normal, projLateral2);

   planeSource->SetOrigin(projOrigin);
   planeSource->SetPoint1(projLateral1);
   planeSource->SetPoint2(projLateral2);

   bool isNormalFlipped = false;
   if (isThick) {
      // At this point, the plane is centered in the middle of the box
      // We create two points to try to guess automatically if we have to flip
      // the sides. We try to see from which border are we closer. That's the tumor
      // border (if we are lucky).
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
      //std::cout << "distFlip: " << distFlip << "; distNoFlip:" << distNotFlip << std::endl;
      double pushDist = 0;
      if (distFlip < distNotFlip)
      {
        isNormalFlipped = true;
        pushDist = -1;
      }
      else
      {
        isNormalFlipped = false;
        pushDist = 1;
      }

      if (forceSideFlip)
      {
          isNormalFlipped = !isNormalFlipped;
          pushDist *= -1;
      }

      planeSource->Push(pushDist);

     // calculates the resolution
     int xresolution = std::sqrt(vtkMath::Distance2BetweenPoints(projOrigin, projLateral1)) / spacing;
     int yresolution = std::sqrt(vtkMath::Distance2BetweenPoints(projOrigin, projLateral2)) / spacing;

     planeSource->SetResolution(xresolution, yresolution);
   }
   planeSource->Update();

   PlaneParameters planeParameters;
   planeParameters.planeSource = planeSource;
   planeParameters.isNormalFlipped = isNormalFlipped;

   return planeParameters;
}
