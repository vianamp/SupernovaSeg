#include <list>
#include <cmath>
#include <cstdio>
#include <random>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <algorithm>

#include <vtkMath.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkLongArray.h>
#include <vtkImageData.h>
#include <vtkImageCast.h>
#include <vtkImageShiftScale.h>
#include <vtkDataSet.h>
#include <vtkDataArray.h>
#include <vtkContourFilter.h>
#include <vtkTransform.h>
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkReverseSense.h>
#include <vtkVectorText.h>
#include <vtkExtractVOI.h>
#include <vtkInformation.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkImageDilateErode3D.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkStructuredPoints.h>
#include <vtkInformationVector.h>
#include <vtkPolyDataNormals.h>
#include <vtkGeometryFilter.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkFillHolesFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkKdTreePointLocator.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkTIFFWriter.h>
#include <vtkPointData.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>
#include <vtkPoints.h>


#include "_database.h"
#include "_Supernova.h"

#define TIFF_ORIENTATION_READER 1

#define DEBUG
