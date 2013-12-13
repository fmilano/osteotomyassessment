cmake_minimum_required(VERSION 2.6)
 
PROJECT(OsteotomyAssessment)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

ADD_EXECUTABLE(OsteotomyAssessment
main.cpp
vtkLORANSACPlane.cxx
vtkOBBTree2.cxx
Projection.cxx)
TARGET_LINK_LIBRARIES(OsteotomyAssessment ${QT_LIBRARIES} vtkCommon vtkIO vtkHybrid)
