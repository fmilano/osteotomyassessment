// Because vtkPlane is strange, we will return as output a PolyData and provide Normal and Origin accessors.
// Author: Federico Milano
// LO-RANSAC (Locally Optimized). Based on RANSAC by David Doria

#ifndef __vtkLORANSACPlane_h
#define __vtkLORANSACPlane_h

#include <vector>

#include <vtkSmartPointer.h>
#include <vtkMinimalStandardRandomSequence.h>

#include <vtkPolyDataAlgorithm.h> //superclass


class vtkPolyData;
class vtkPlane;
class vtkPoints;

class vtkInformation;
class vtkInformationVector;

class vtkLORANSACPlane : public vtkPolyDataAlgorithm
{
public:
   static vtkLORANSACPlane *New();
   vtkTypeRevisionMacro(vtkLORANSACPlane, vtkPolyDataAlgorithm);
   void PrintSelf(ostream &os, vtkIndent indent);

   std::vector<unsigned int> MaxInlierIndices;

   // Description:
   // Specify the source object. This is the object that will be moved during the transformation.
   vtkPolyData *GetSource();

   // Description:
   // Specify the target object. This is the object that will stay in place.
   vtkPolyData *GetTarget();

   void AddSourceConnection(vtkAlgorithmOutput* input);
   void RemoveAllSources();

   vtkSmartPointer<vtkPlane> GetBestPlane() { return BestPlane; }

   vtkPolyData* GetOutput();

   void SetBounds(double bounds[]);

  //double RandomPoint0[3];
  //double RandomPoint1[3];
  //double RandomPoint2[3];

  vtkSetMacro(InlierThreshold, double);
  vtkGetMacro(InlierThreshold, double);
      
  vtkSetMacro(MaxIterations, unsigned int);
  vtkGetMacro(MaxIterations, unsigned int);
        
  vtkSetMacro(GoodEnough, double);
  vtkGetMacro(GoodEnough, double);
  
 protected:
   vtkLORANSACPlane();
   ~vtkLORANSACPlane();

  int FillInputPortInformation( int port, vtkInformation* info );
    
  // Generate output
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

 private:
   double InlierThreshold; //how close to the plane a point must be to be considered an inlier
   unsigned int MaxIterations; //how many times to pick points and fit a plane
   unsigned int NumPointsToFit; //how many points to pick and fit a plane to in each iteration
   double GoodEnough; //The percentage of the data that is fit to consider the model "good enough"
   
   std::vector<unsigned int> DetermineInliers(vtkPoints* points, vtkPlane* plane);
   std::vector<unsigned int> UniqueRandomIndices(const unsigned int maxIndex, const unsigned int numIndices);
   
   vtkSmartPointer<vtkPlane> BestPlane;

   vtkSmartPointer<vtkMinimalStandardRandomSequence> randomSequence;

   double Bounds[6];
};

////////// Helper Functions /////////////
void ExtractPoints(vtkPointSet* polydata, const std::vector<unsigned int>& indices, vtkPoints* output);
void BestFitPlane(vtkPoints *points, vtkPlane *BestPlane);
void CopyPlane(vtkPlane* plane, vtkPlane* output);
void CenterOfMass(vtkPoints* points, double* center);
    
#endif
