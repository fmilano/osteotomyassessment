#include <algorithm>
#include <cmath>

#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
//#include "vtkCommand.h"

#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkPlane.h"
#include "vtkSmartPointer.h"
#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkPlaneSource.h"
#include "vtkTriangle.h"

#include <set>

#include "vtkLORANSACPlane.h"

vtkStandardNewMacro(vtkLORANSACPlane);
vtkCxxRevisionMacro(vtkLORANSACPlane, "$Revision$");


//-----------------------------------------------------------------------------
vtkLORANSACPlane::vtkLORANSACPlane()
{
  this->InlierThreshold = 1.0;
  this->MaxIterations = 1000;
  this->NumPointsToFit = 3; 
  this->GoodEnough = 1.0; 
  
  this->SetNumberOfInputPorts( 1 );
  this->SetNumberOfOutputPorts( 1 );
}


//-----------------------------------------------------------------------------
vtkLORANSACPlane::~vtkLORANSACPlane()
{
}


//----------------------------------------------------------------------------
void vtkLORANSACPlane::AddSourceConnection(vtkAlgorithmOutput* input)
{
  this->AddInputConnection(1, input);
}


//----------------------------------------------------------------------------
void vtkLORANSACPlane::RemoveAllSources()
{
  this->SetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
vtkPolyData* vtkLORANSACPlane::GetOutput()
{
  return vtkPolyData::SafeDownCast(this->GetOutputDataObject(0));
}

//----------------------------------------------------------------------------
int vtkLORANSACPlane::FillInputPortInformation( int port, vtkInformation* info )
{
  if (!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }
  if ( port == 0 )
    {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet" );
    return 1;
    }
  return 0;
}

void vtkLORANSACPlane::SetBounds(double bounds[]) {
   Bounds[0] = bounds[0];
   Bounds[1] = bounds[1];
   Bounds[2] = bounds[2];
   Bounds[3] = bounds[3];
   Bounds[4] = bounds[4];
   Bounds[5] = bounds[5];
}

//----------------------------------------------------------------------------
int vtkLORANSACPlane::RequestData(vtkInformation *vtkNotUsed(request),
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkPointSet *input = vtkPointSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    //track best model
    //vtkSmartPointer<vtkPlane> bestPlane = vtkSmartPointer<vtkPlane>::New();
    BestPlane = vtkSmartPointer<vtkPlane>::New();

    //track number of inliers of best model
    unsigned int maxInliers = 0;

    //seed the random number gererator
    // !!! ///

    for (unsigned int iter = 0; iter < MaxIterations; iter++) {
        //pick NumPointsToFit random indices
        std::vector<unsigned int> randomIndices = UniqueRandomIndices(input->GetNumberOfPoints(), NumPointsToFit);

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        ExtractPoints(input, randomIndices, points);

        //find the best plane through these random points
        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        //BestFitPlane(points, plane);
        // -------------------------------
        // Calculates the best fit plane
        double p0[3], p1[3], p2[3];
        input->GetPoint(randomIndices[0], p0);
        input->GetPoint(randomIndices[1], p1);
        input->GetPoint(randomIndices[2], p2);

        double normal[3];
        vtkTriangle::ComputeNormal(p0, p1, p2, normal);

        plane->SetNormal(normal);
        plane->SetOrigin(p0);

        // -------------------------------

        std::vector<unsigned int> inlierIndices = DetermineInliers(input->GetPoints(), plane);

        if (inlierIndices.size() > maxInliers) {

            // save model
            maxInliers = inlierIndices.size();
            CopyPlane(plane, BestPlane);
            MaxInlierIndices = inlierIndices;

            // Run local optimization
            vtkSmartPointer<vtkPoints> inliers = vtkSmartPointer<vtkPoints>::New();
            ExtractPoints(input, inlierIndices, inliers);

            BestFitPlane(inliers, plane);
            inlierIndices = DetermineInliers(input->GetPoints(), plane);

            // save optimized model
            if (inlierIndices.size() > maxInliers) {
                std::cerr << "Optimizing!" << std::endl;

                maxInliers = inlierIndices.size();
                CopyPlane(plane, BestPlane);
                MaxInlierIndices = inlierIndices;
            }
        }

        if(inlierIndices.size() > input->GetNumberOfPoints() * GoodEnough) //if GoodEnough % of the points fit the model, we can stop the search
        {
            break;
        }

    } //end ransac loop

    /*cout << "Best plane: " << endl << "---------------" << endl;
    cout << "Total points: " << input->GetNumberOfPoints() << endl;
    cout << "Inliers: " << maxInliers << endl;
    */
    double n[3];
    BestPlane->GetNormal(n);

    //double bounds[6];
    //input->ComputeBounds();
    //input->GetBounds(bounds);
    double minBound[3];
    minBound[0] = Bounds[0];
    minBound[1] = Bounds[2];
    minBound[2] = Bounds[4];
    double maxBound[3];
    maxBound[0] = Bounds[1];
    maxBound[1] = Bounds[3];
    maxBound[2] = Bounds[5];

    int index = 0;
    if (fabs(Bounds[index] - Bounds[index + 1]) < fabs(Bounds[2] - Bounds[3])) {
        index = 2;
    }
    if (fabs(Bounds[index] - Bounds[index + 1]) < fabs(Bounds[4] - Bounds[5])) {
        index = 4;
    }

    double origin[3];
    origin[0] = Bounds[1];
    origin[1] = Bounds[3];
    origin[2] = Bounds[5];

    origin[index / 2] = Bounds[index];

    double projOrigin[3], projMin[3], projMax[3];

    BestPlane->ProjectPoint(origin, projOrigin);
    BestPlane->ProjectPoint(minBound, projMin);
    BestPlane->ProjectPoint(maxBound, projMax);

    vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
    planeSource->SetNormal(n);
    planeSource->SetOrigin(projOrigin);
    planeSource->SetPoint1(projMax);
    planeSource->SetPoint2(projMin);

    int xresolution = std::sqrt(vtkMath::Distance2BetweenPoints(projOrigin, projMax)) / 0.5;
    int yresolution = std::sqrt(vtkMath::Distance2BetweenPoints(projOrigin, projMin)) / 0.5;

    planeSource->SetResolution(xresolution, yresolution);
    planeSource->Update();

    output->DeepCopy(planeSource->GetOutput());

    return 1;
}

void vtkLORANSACPlane::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

std::vector<unsigned int> vtkLORANSACPlane::DetermineInliers(vtkPoints* points, vtkPlane* plane)
{
    //find the distance from every point to the plane
    std::vector<unsigned int> inlierIndices;
    for(unsigned int i = 0; i < points->GetNumberOfPoints(); i++) {
        //double distance = fabs(vgl_distance(Plane, Points[i]));
        double point[3];
        points->GetPoint(i, point);
        //double distance = plane->DistanceToPlane(p);
        double n[3];
        double o[3];
        plane->GetNormal(n);
        plane->GetOrigin(o);

        double distance = vtkPlane::DistanceToPlane(point, n, o);

        if(distance < this->InlierThreshold)
        {
           inlierIndices.push_back(i);
        }
    }

    return inlierIndices;
}
    
///////////////////////////////////////
////////// Helper Functions /////////////
////////////////////////////////////////
std::vector<unsigned int> UniqueRandomIndices(const unsigned int maxIndex, const unsigned int numIndices)
{
  //generate Number unique random indices from 0 to MAX

   std::vector<unsigned int> indices;
  
   //SeedRandom();
   //cannot generate more unique numbers than than the size of the set we are sampling
   if(!(numIndices <= maxIndex+1)) {
      return indices;
   }
		
	std::set<unsigned int> S;
   while(S.size() < numIndices)
   {
      S.insert(vtkMath::Random(0, maxIndex));
   }

   for(std::set<unsigned int>::iterator iter = S.begin(); iter != S.end(); iter++)
   {
      indices.push_back(*iter);
   }

   return indices;
}

void ExtractPoints(vtkPointSet* points, const std::vector<unsigned int>& indices, vtkPoints* output)
{
    for(unsigned int i = 0; i < indices.size(); i++) {
        double p[3];
        points->GetPoint(indices[i], p);
        output->InsertNextPoint(p[0], p[1], p[2]);
    }
}


/* allocate memory for an nrow x ncol matrix */
template<class TReal>
    TReal **create_matrix ( long nrow, long ncol )
{
  typedef TReal* TRealPointer;
  TReal **m = new TRealPointer[nrow];

  TReal* block = ( TReal* ) calloc ( nrow*ncol, sizeof ( TReal ) );
  m[0] = block;
  for ( int row = 1; row < nrow; ++row )
  {
    m[ row ] = &block[ row * ncol ];
  }
  return m;
}

/* free a TReal matrix allocated with matrix() */
template<class TReal>
    void free_matrix ( TReal **m )
{
  free ( m[0] );
  delete[] m;
}

void BestFitPlane(vtkPoints *points, vtkPlane *BestPlane)
{
  //Compute the best fit (least squares) plane through a set of points.
  vtkIdType NumPoints = points->GetNumberOfPoints();
  double dNumPoints = static_cast<double>(NumPoints);
  
  //find the center of mass of the points
  double Center[3];
  CenterOfMass(points, Center);
  //vtkstd::cout << "Center of mass: " << Center[0] << " " << Center[1] << " " << Center[2] << vtkstd::endl;
  
  //Compute sample covariance matrix
  double **a = create_matrix<double> ( 3,3 );
  a[0][0] = 0; a[0][1] = 0;  a[0][2] = 0;
  a[1][0] = 0; a[1][1] = 0;  a[1][2] = 0;
  a[2][0] = 0; a[2][1] = 0;  a[2][2] = 0;
  
  for(unsigned int pointId = 0; pointId < NumPoints; pointId++ )
  {
    double x[3];
    double xp[3];
    points->GetPoint(pointId, x);
    xp[0] = x[0] - Center[0]; 
    xp[1] = x[1] - Center[1]; 
    xp[2] = x[2] - Center[2];
    for (unsigned int i = 0; i < 3; i++)
    {
      a[0][i] += xp[0] * xp[i];
      a[1][i] += xp[1] * xp[i];
      a[2][i] += xp[2] * xp[i];
    }
  }
  
    //divide by N-1
  for(unsigned int i = 0; i < 3; i++)
  {
    a[0][i] /= dNumPoints-1;
    a[1][i] /= dNumPoints-1;
    a[2][i] /= dNumPoints-1;
  }

  // Extract eigenvectors from covariance matrix
  double **eigvec = create_matrix<double> ( 3,3 );
  
  double eigval[3];
  vtkMath::Jacobi(a,eigval,eigvec);

    //Jacobi iteration for the solution of eigenvectors/eigenvalues of a 3x3 real symmetric matrix. Square 3x3 matrix a; output eigenvalues in w; and output eigenvectors in v. Resulting eigenvalues/vectors are sorted in decreasing order; eigenvectors are normalized.
  
  //Set the plane normal to the smallest eigen vector
  BestPlane->SetNormal(eigvec[0][2], eigvec[1][2], eigvec[2][2]);
  
  //cleanup
  free_matrix(eigvec);
  free_matrix(a);
  
  //Set the plane origin to the center of mass
  BestPlane->SetOrigin(Center[0], Center[1], Center[2]);
}

void CopyPlane(vtkPlane* plane, vtkPlane* output)
{
  double n[3];
  plane->GetNormal(n);
  
  double o[3];
  plane->GetOrigin(o);
  
  output->SetNormal(n);
  output->SetOrigin(o);
  
}

void CenterOfMass(vtkPoints* points, double* center)
{
  center[0] = 0.0;
  center[1] = 0.0;
  center[2] = 0.0;
    
  for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
  {
    double point[3];
    points->GetPoint(i, point);
    
    center[0] += point[0];
    center[1] += point[1];
    center[2] += point[2];
  }
  
  double numberOfPoints = static_cast<double>(points->GetNumberOfPoints());
  center[0] = center[0]/numberOfPoints;
  center[1] = center[1]/numberOfPoints;
  center[2] = center[2]/numberOfPoints;
}


