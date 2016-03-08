#ifndef _SUPERNOVA_H
#define _SUPERNOVA_H

#include "includes.h"

    struct _point {
        double x, y, z;
    };

    class _Supernova {
        
        private:    

            int _rmax, _nrays, _freq;
            vtkPolyData *Rays, *Peaks, *Cell;
            double _xo, _yo, _zo, _scalefactor;
        
        public:

            int GetMaximumRadius() {
                return _rmax;
            }    
            void SetMaximumRadius(int rmax) {
                _rmax = rmax;
            }
            void SetCenter(double xo, double yo, double zo) {
                _xo = xo;
                _yo = yo;
                _zo = zo;
            }
            void SetScaleFactor(double _dxy, double _dz) {
                _scalefactor = _dz/_dxy;
            }
            int GetNumberOfRays() {
                return _nrays;
            }
            void SetNumberOfRays(int nrays) {
                _nrays = nrays;
            }
            void SetFrequency(int freq) {
                _freq = freq;
            }
            void SetPoints(vtkPoints *Points) {
                Rays -> SetPoints(Points);
            }
            void SetLines(vtkCellArray *Array) {
                Rays -> SetLines(Array);
            }
            vtkDataArray *GetIntensities() {
                return Rays -> GetPointData() -> GetScalars();
            } 
            void SetLinesIntensity(vtkUnsignedShortArray *Intensities) {
                Rays -> GetPointData() -> SetScalars(Intensities);
            }
            void Initialize();

            void Probe(vtkImageData *Image);

            void Save(const char FileName[], const int _id);

            void ApplyLimits(const double r1, const bool force);

            void Segmentation();

            void ClipImageData(const char MitoFileName[], vtkImageData *Image, vtkImageData *ClipImage, const int _id);

            void ScalePolyData(double _dxy, double _dz);

            void AjustCoordinates(int Ly, double _dxy, double _dz);

            void EstimateImageMeanAndStd(vtkImageData *Image, double *_mean, double *_std);

            void GetXYZFromRay(const int ray, double *x, double *y, double *z);

            _Supernova() {
                _rmax = 55;
                _freq = 20;
                _nrays = 1000;
            }

            ~_Supernova();

    };

#endif