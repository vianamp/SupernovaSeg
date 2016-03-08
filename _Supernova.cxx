#include "_Supernova.h"

    int ssdx_sort[26] = { 0,-1, 0, 1, 0, 0,-1, 0, 1, 0,-1, 1, 1,-1,-1, 0, 1, 0, -1, 1, 1,-1,-1, 1, 1,-1};
    int ssdy_sort[26] = { 0, 0,-1, 0, 1, 0, 0,-1, 0, 1,-1,-1, 1, 1, 0,-1, 0, 1, -1,-1, 1, 1,-1,-1, 1, 1};
    int ssdz_sort[26] = {-1, 0, 0, 0, 0, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1, -1,-1,-1,-1, 1, 1, 1, 1};

    void _Supernova::GetXYZFromRay(const int ray, double *x, double *y, double *z) {
        double t = 2.0 * 3.1415*(1.0*ray)/_nrays;
        *z = 3.1415 - t;
        *x = cos(_freq*t);
        *y = sin(_freq*t);
    }

    void _Supernova::Initialize() {

        int m, p;
        double i, j, k;
        double t, x, y, z, n;
        vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> RayArray = vtkSmartPointer<vtkCellArray>::New();

        vtkIdType id = 0;
        for (m = 0; m < _nrays; m++) {

            GetXYZFromRay(m,&x,&y,&z);

            n = sqrt(x*x+y*y+z*z);

            RayArray -> InsertNextCell(_rmax);
            for (p = 0; p < _rmax; p++) {
                i = _xo + p * (x / n);
                j = _yo + p * (y / n);
                k = _zo + p * (z / n);
                Points -> InsertNextPoint(i,j,k);
                RayArray -> InsertCellPoint(id);
                id++;
            }

        }

        Rays = vtkPolyData::New();
        Rays -> SetPoints(Points);
        Rays -> SetLines(RayArray);
        
    }

    void _Supernova::Probe(vtkImageData *Image) {
        double v;
        int i, j, k;
        vtkIdType id;
        int *Dim = Image -> GetDimensions();

        vtkSmartPointer<vtkCellArray> RayArray = vtkSmartPointer<vtkCellArray>::New();

        vtkSmartPointer<vtkUnsignedShortArray> Intensities = vtkSmartPointer<vtkUnsignedShortArray>::New();
        Intensities -> SetNumberOfComponents(1);
        Intensities -> SetNumberOfTuples(_nrays*_rmax);
        Intensities -> FillComponent(0,0);

        double r[3];
        for (id = 0; id < Rays->GetNumberOfPoints(); id++) {
            v = 0;
            Rays -> GetPoint(id,r);
            i = (int)r[0]; j = (int)r[1]; k = (int)r[2];
            if ( (i >= 0 && i < Dim[0]) && (j >= 0 && j < Dim[1]) && (k >= 0 && k < Dim[2]) ) {
                v = Image -> GetScalarComponentAsDouble(i,j,k,0);
            }
            if (k<5) v = 0; // THIS LINE IS CRITICAL FOR GOOD CELL SHAPE IN Z
                            // MIGHT HAVE TO DO A AUTOMATIC DETECTION OF 5
            Intensities -> SetTuple1(id,v);
        }

        Rays -> GetPointData() -> SetScalars(Intensities);

        #ifdef DEBUG
            vtkSmartPointer<vtkPolyDataWriter> Writer =  vtkSmartPointer<vtkPolyDataWriter>::New();
            Writer -> SetFileTypeToBinary();
            Writer -> SetInputData(Rays);
            Writer -> SetFileName("Rays.vtk");
            Writer -> Write();
            printf("Rays saved in Supernove folder under the name Rays.vtk!\n");
        #endif

    }

    void _Supernova::ApplyLimits(const double r1, bool force) {
        vtkIdType r, i, id;
        double v, x, y, z, w;
        vtkSmartPointer<vtkDataArray> Scalars = Rays -> GetPointData() -> GetScalars();
        for (r = 0; r < _nrays; r++) {
            GetXYZFromRay(r,&x,&y,&z);
            w = r1 / ( fabs(z/sqrt(x*x+y*y+z*z)) * (_scalefactor-1) + 1 );
            for (i = 0; i < _rmax; i++) {
                id = i + r * _rmax;
                if (force) {
                    if (i == (int)w) Scalars -> SetTuple1(id,65535);
                } else {
                    v = Scalars -> GetTuple1(id);
                    v = (i<w)?v:v/0.0;
                    Scalars -> SetTuple1(id,v);
                }
            }
        }
    }

    void _Supernova::Save(const char FileName[], const int _id) {

        vtkSmartPointer<vtkPolyDataWriter> SurfaceWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
        SurfaceWriter -> SetFileTypeToBinary();
        SurfaceWriter -> SetInputData(Cell);
        SurfaceWriter -> SetFileName(FileName);
        SurfaceWriter -> Write();

    }

    void _Supernova::Segmentation() {

       #ifdef DEBUG
            printf("Saving Max projection...\n");
        #endif

        vtkSmartPointer<vtkImageData> Proj = vtkSmartPointer<vtkImageData>::New();
        Proj -> SetDimensions(_rmax,_nrays,1);
        
        vtkSmartPointer<vtkDataArray> Intensities = Rays -> GetPointData() -> GetScalars();
        
        #ifdef DEBUG
            Proj -> GetPointData() -> SetScalars(Intensities);
            vtkSmartPointer<vtkTIFFWriter> TIFFWriter = vtkSmartPointer<vtkTIFFWriter>::New();
            TIFFWriter -> SetFileName("Projection.tif");
            TIFFWriter -> SetFileDimensionality(2);
            TIFFWriter -> SetInputData(Proj);
            TIFFWriter -> Write();
            printf("Projection saved in Supernove folder under the name Projection.tif!\n");
        #endif

        int i, j, m;
        std::vector<double> Prec;
        std::vector<double>::iterator it;

        vtkSmartPointer<vtkDoubleArray> Map = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkUnsignedShortArray> Step = vtkSmartPointer<vtkUnsignedShortArray>::New();
        Map -> DeepCopy(Intensities);
        Step -> DeepCopy(Intensities);
        Step -> FillComponent(0,0);

        for (i = 1; i < _nrays; i++) {
            for (j = 0; j < _rmax; j++) {
                for (m = -1; m <=1; m++) {
                    if ( (j+m >= 0) && (j+m < _rmax) ) {
                        Prec.push_back(Map->GetTuple1((j+m)+(i-1)*_rmax));
                    }
                }
                it = std::max_element(Prec.begin(),Prec.end());
                Step -> SetTuple1(j+i*_rmax,(int)std::distance(Prec.begin(),it));
                Map -> SetTuple1(j+i*_rmax,Intensities->GetTuple1(j+i*_rmax)+(*it));
                Prec.clear();
            }
        }

        Map -> Modified();
        Step -> Modified();

        double v, map_range[2];
        Map -> GetRange(map_range);
        printf("Range = %1.3f / %1.3f\n",map_range[0],map_range[1]);

        for (i = 0; i < Intensities->GetNumberOfTuples(); i++) {
            v = 65535*(Map->GetTuple1(i)-map_range[0])/(map_range[1]-map_range[0]);
            Intensities -> SetTuple1(i,v);
            Map -> SetTuple1(i,v);
        }

        Intensities -> Modified();

        #ifdef DEBUG
            Proj -> GetPointData() -> SetScalars(Intensities);
            TIFFWriter = vtkSmartPointer<vtkTIFFWriter>::New();
            TIFFWriter -> SetFileName("Map.tif");
            TIFFWriter -> SetFileDimensionality(2);
            TIFFWriter -> SetInputData(Proj);
            TIFFWriter -> Write();
            printf("Map saved in Supernove folder under the name Map.tif!\n");
        #endif

        #ifdef DEBUG
            Proj -> GetPointData() -> SetScalars(Step);
            TIFFWriter = vtkSmartPointer<vtkTIFFWriter>::New();
            TIFFWriter -> SetFileName("Step.tif");
            TIFFWriter -> SetFileDimensionality(2);
            TIFFWriter -> SetInputData(Proj);
            TIFFWriter -> Write();
            printf("Steps saved in Supernove folder under the name Step.tif!\n");
        #endif

        int ext = 0;
        long int jo = 0;
        i = _nrays - 1;
        std::vector< std::pair<int,int> > Path;
        for (j = 0; j < _rmax; j++) {
            if (Intensities->GetTuple1(j+i*_rmax) == 65535) {
                ext++;
                jo += j;
            }
        }
        int delta;
        jo /= ext;
        Intensities -> FillComponent(0,0);
        do {
            Path.push_back(std::make_pair(i,jo));
            Intensities -> SetTuple1(jo+i*_rmax,65535);
            delta = Step -> GetTuple1(jo+i*_rmax) - 1;
            jo += delta;//(jo+delta>=0) ? delta : 0;
            i--;
        } while (i >= 0);


        #ifdef DEBUG
            Proj -> GetPointData() -> SetScalars(Intensities);
            TIFFWriter = vtkSmartPointer<vtkTIFFWriter>::New();
            TIFFWriter -> SetFileName("Path.tif");
            TIFFWriter -> SetFileDimensionality(2);
            TIFFWriter -> SetInputData(Proj);
            TIFFWriter -> Write();
            printf("Path saved in Supernove folder under the name Path.tif!\n");
        #endif

        vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();

        int p, a, b;
        a = Path.front().second;
        b = Path.back().second;
        double t, x, y, z, n;
        do {
            m = Path.back().first;
            p = Path.back().second;
            GetXYZFromRay(m,&x,&y,&z);            
            n = sqrt(x*x+y*y+z*z);
            Points -> InsertNextPoint(_xo + p * (x / n),_yo + p * (y / n),_zo + p * (z / n));
            Path.pop_back();
        } while (Path.size() > 0);

        Peaks = vtkPolyData::New();
        Peaks -> SetPoints(Points);

        vtkSmartPointer<vtkPolyDataWriter> Writer;

        #ifdef DEBUG
            Writer =  vtkSmartPointer<vtkPolyDataWriter>::New();
            Writer -> SetFileTypeToBinary();
            Writer -> SetInputData(Peaks);
            Writer -> SetFileName("Peaks.vtk");
            Writer -> Write();
            printf("Peaks saved in Supernove folder under the name Peaks.vtk!\n");
        #endif

        vtkSmartPointer<vtkSurfaceReconstructionFilter> PoissonRecFilter = vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
        PoissonRecFilter -> SetNeighborhoodSize(200);
        PoissonRecFilter -> SetInputData(Peaks);
        PoissonRecFilter -> Update();

        vtkSmartPointer<vtkContourFilter> ContourFilter = vtkSmartPointer<vtkContourFilter>::New();
        ContourFilter -> SetInputConnection(PoissonRecFilter->GetOutputPort());
        ContourFilter -> SetValue(0, 0.0);
        ContourFilter -> Update();

        Cell = vtkPolyData::New();
        Cell -> DeepCopy(ContourFilter->GetOutput());

        #ifdef DEBUG
            Writer =  vtkSmartPointer<vtkPolyDataWriter>::New();
            Writer -> SetFileTypeToBinary();
            Writer -> SetInputData(Cell);
            Writer -> SetFileName("Cell.vtk");
            Writer -> Write();
            printf("Cell saved in Supernove folder under the name Cell.vtk!\n");
        #endif

    }

    void FillHoles(vtkSmartPointer<vtkImageData> ImageData) {

        #ifdef DEBUG
            printf("\tSearching for holes in the image...\n");
        #endif

        int *Dim = ImageData -> GetDimensions();
        vtkIdType N = ImageData -> GetNumberOfPoints();

        vtkIdType i, s, ido, id;

        int x, y, z;
        double v, r[3];
        bool find = true;
        long long int ro[3];
        long int scluster, label;
        ro[0] = Dim[0] * Dim[1] * Dim[2];
        ro[1] = Dim[0] * Dim[1];

        vtkSmartPointer<vtkIdList> CurrA = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> NextA = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkLongArray> CSz = vtkSmartPointer<vtkLongArray>::New();
        vtkSmartPointer<vtkLongArray> Volume = vtkSmartPointer<vtkLongArray>::New();
        Volume -> SetNumberOfComponents(1);
        Volume -> SetNumberOfTuples(N);
        Volume -> FillComponent(0,0);

        for (x = 1; x < Dim[0]-1; x++) {
            for (y = 1; y < Dim[1]-1; y++) {
                for (z = 1; z < Dim[2]-1; z++) {
                    id = ImageData -> FindPoint(x,y,z);
                    if ((unsigned short int)ImageData->GetScalarComponentAsDouble(x,y,z,0)) {
                        Volume -> SetTuple1(id,1);
                    } else {
                        Volume -> SetTuple1(id,0);
                    }
                }
            }
        }

        Volume -> Modified();

        label = 0;
        while (find) {
            for (s = 0; s < CurrA->GetNumberOfIds(); s++) {
                ido = CurrA -> GetId(s);
                ImageData -> GetPoint(ido,r);
                x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
                for (i = 0; i < 6; i++) {
                    id = ImageData -> FindPoint(x+ssdx_sort[i],y+ssdy_sort[i],z+ssdz_sort[i]);
                    v = Volume -> GetTuple1(id);
                    if ((long int)v > 0) {
                        NextA -> InsertNextId(id);
                        Volume -> SetTuple1(id,-label);
                        scluster++;
                    }
                }
            }
            if (!NextA->GetNumberOfIds()) {
                find = false;
                for (id=ro[0]; id--;) {
                    v = Volume -> GetTuple1(id);
                    if ((long int)v > 0) {
                        find = true;
                        ro[0] = id;
                        break;
                    }
                }
                if (label) {
                    CSz -> InsertNextTuple1(scluster);
                }
                if (find) {
                    label++;
                    scluster = 1;
                    Volume -> SetTuple1(id,-label);
                    CurrA -> InsertNextId(id);
                }
            } else {
                CurrA -> Reset();
                CurrA -> DeepCopy(NextA);
                NextA -> Reset();
            }
        }

        int biggest;
        long int max_size = 0;
        for (id = 0; id < CSz->GetNumberOfTuples(); id++) {
            if (CSz->GetTuple1(id) > max_size) {
                biggest = id;
                max_size = CSz->GetTuple1(id);
            }
        }

        biggest  = -(biggest+1);

        printf("Biggest = %d\n",biggest);

        for (id = N; id--;) {
            if ((long int)Volume->GetTuple1(id) == biggest) {
                ImageData -> GetPointData() -> GetScalars() -> SetTuple1(id,255);
            } else if ((long int)Volume->GetTuple1(id) < biggest) {
                ImageData -> GetPointData() -> GetScalars() -> SetTuple1(id,0);
            }
        }
        ImageData -> GetPointData() -> GetScalars() -> Modified();

        #ifdef DEBUG
            printf("\tNumber of filled holes: %ld\n",(long int)CSz->GetNumberOfTuples()-1);
        #endif

    }

    void _Supernova::EstimateImageMeanAndStd(vtkImageData *Image, double *_mean, double *_std) {
        std::vector<unsigned short> v;
        for (vtkIdType id = Image->GetNumberOfPoints(); id--;) {
            v.push_back(Image->GetPointData()->GetScalars()->GetTuple1(id));
        }
        double sum = std::accumulate(v.begin(), v.end(), 0.0);
        double mean = sum / v.size();
        std::vector<double> diff(v.size());
        std::transform(v.begin(), v.end(), diff.begin(),std::bind2nd(std::minus<double>(), mean));
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / v.size());
        *_mean = mean;
        *_std = stdev;
        printf("\t <I> = %1.3f\n",mean);
        printf("\t sd_I = %1.3f\n",stdev);
    }

    void _Supernova::ClipImageData(const char MitoFileName[], vtkImageData *Image, vtkImageData *ClipImage, const int _id) {

        #ifdef DEBUG
            printf("Clipping 3D Volume...\n");
        #endif

        vtkSmartPointer<vtkSelectEnclosedPoints> Filter = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
        Filter -> SetInputData(Image);
        Filter -> SetSurfaceData(Cell);
        Filter -> Update();

        double range[2];
        Filter -> GetOutput() -> GetScalarRange(range);

        printf("Range = [%1.3f,%1.3f]\n",range[0],range[1]);

        vtkSmartPointer<vtkImageShiftScale> Cast = vtkSmartPointer<vtkImageShiftScale>::New();
        Cast -> SetInputData(Filter->GetOutput());
        Cast -> SetShift(0.0);
        Cast -> SetScale(65535.0/range[1]);
        Cast -> SetOutputScalarTypeToUnsignedShort();
        Cast -> Update();

        vtkSmartPointer<vtkImageData> BinaryImage = Cast -> GetOutput();

        unsigned short *qPixel, *pPixel = static_cast<unsigned short*>(BinaryImage->GetScalarPointer());
        for( vtkIdType id = BinaryImage->GetNumberOfPoints(); id--; )
            if (pPixel[id]>0) pPixel[id] = 65535;

        vtkSmartPointer<vtkImageGaussianSmooth> Gauss = vtkSmartPointer<vtkImageGaussianSmooth>::New();
        Gauss -> SetInputData(BinaryImage);
        Gauss -> SetStandardDeviations(3,3,2);
        Gauss -> SetRadiusFactors(3,3,2);
        Gauss -> Update();

        vtkSmartPointer<vtkImageData> SmoothImage = vtkSmartPointer<vtkImageData>::New();
        SmoothImage -> DeepCopy(Gauss->GetOutput());

        int i, j, *Dim = BinaryImage -> GetDimensions();
        for (i = 0; i < Dim[0]; i++) {
            for (j = 0; j < Dim[1]; j++) {
                SmoothImage -> SetScalarComponentFromDouble(i,j,0,0,0);
                SmoothImage -> SetScalarComponentFromDouble(i,j,Dim[2]-1,0,0);
            }
        }

        double threshold = 0.5*65535;
        vtkSmartPointer<vtkContourFilter> ContourFilter = vtkSmartPointer<vtkContourFilter>::New();
        ContourFilter -> SetInputData(SmoothImage);
        ContourFilter -> SetValue(0,threshold);
        ContourFilter -> Update();

        Cell -> DeepCopy(ContourFilter->GetOutput());

        vtkSmartPointer<vtkIntArray> ID = vtkSmartPointer<vtkIntArray>::New();
        ID -> SetName("ID");
        ID -> SetNumberOfTuples(Cell->GetNumberOfPoints());
        ID -> SetNumberOfComponents(1);
        ID -> FillComponent(0,_id);

        Cell -> GetPointData() -> SetScalars(ID);
        Cell -> Modified();

        vtkSmartPointer<vtkTIFFReader> Reader = vtkSmartPointer<vtkTIFFReader>::New();
        Reader -> SetFileName(MitoFileName);
        Reader -> SetOrientationType(TIFF_ORIENTATION_READER);
        Reader -> Update();

        vtkSmartPointer<vtkImageData> Mito = Reader -> GetOutput();

        double mean, stdev;
        EstimateImageMeanAndStd(Mito,&mean,&stdev);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> NormalDist(mean,0.5*stdev);

        int voi[6];
        double r[3];
        voi[0] = voi[2] = voi[4] = 1E3;
        voi[1] = voi[3] = voi[5] = 0;
        qPixel = static_cast<unsigned short*>(Mito->GetScalarPointer());
        pPixel = static_cast<unsigned short*>(SmoothImage->GetScalarPointer());
        for( vtkIdType id = SmoothImage->GetNumberOfPoints(); id--; ) {
            if (pPixel[id]>threshold) {
                BinaryImage -> GetPoint(id,r);
                voi[0] = ((int)r[0] < voi[0]) ? (int)r[0] : voi[0];
                voi[1] = ((int)r[0] > voi[1]) ? (int)r[0] : voi[1];
                voi[2] = ((int)r[1] < voi[2]) ? (int)r[1] : voi[2];
                voi[3] = ((int)r[1] > voi[3]) ? (int)r[1] : voi[3];
                voi[4] = ((int)r[2] < voi[4]) ? (int)r[2] : voi[4];
                voi[5] = ((int)r[2] > voi[5]) ? (int)r[2] : voi[5];
            } else {
                qPixel[id] = (unsigned short)NormalDist(gen);
            }
        }

        voi[0] -= (voi[0]-5>0) ? 5 : 0;
        voi[2] -= (voi[2]-5>0) ? 5 : 0;
        voi[1] += (voi[1]+5<Dim[0]) ? 5 : 0;
        voi[3] += (voi[3]+5<Dim[1]) ? 5 : 0;

        printf("VOI: x = [%d,%d], y = [%d,%d], z = [%d,%d]\n",voi[0],voi[1],voi[2],voi[3],0,Dim[2]-1);

        vtkSmartPointer<vtkExtractVOI> ExtractVOI = vtkSmartPointer<vtkExtractVOI>::New();
        ExtractVOI -> SetVOI(voi[0],voi[1],voi[2],voi[3],voi[4],voi[5]);
        ExtractVOI -> SetSampleRate(1,1,1);
        ExtractVOI -> SetInputData(Mito);
        ExtractVOI -> Update();

        ClipImage -> DeepCopy(ExtractVOI->GetOutput());

        #ifdef DEBUG
            printf("\t>Done.\n");
        #endif

    }

    void _Supernova::AjustCoordinates(int Ly, double _dxy, double _dz) {
        vtkSmartPointer<vtkCleanPolyData> Clean = vtkSmartPointer<vtkCleanPolyData>::New();
        Clean -> SetInputData(Cell);
        Clean -> Update();

        Cell -> DeepCopy(Clean -> GetOutput());

        vtkIdType i;
        double r[3], x, y, z;
        double xcm, ycm, zcm;
        xcm = ycm = zcm = 0.0;
        for (i = 0; i < Cell -> GetNumberOfPoints(); i++) {
            Cell -> GetPoint(i,r);
            xcm += r[0];
            ycm += r[1];
            zcm += r[2];
        }
        xcm /= Cell -> GetNumberOfPoints();
        ycm /= Cell -> GetNumberOfPoints();
        zcm /= Cell -> GetNumberOfPoints();

        printf("%1.3f\t%1.3f\t%1.3f\n",xcm,ycm,zcm);

        for (i = 0; i < Cell -> GetNumberOfPoints(); i++) {
            Cell -> GetPoint(i,r);
            x = 1.4 * _dxy * ( r[0] - xcm ) + _dxy * xcm;
            y = 1.4 * _dxy * ( r[1] - ycm ) + _dxy * ycm;
            z = 1.4 * _dz  * ( r[2] - zcm ) + _dz  * zcm;
            Cell -> GetPoints() -> SetPoint(i,x,_dxy*Ly-y,z);
        }
        Cell -> Modified();

        vtkSmartPointer<vtkFillHolesFilter> Holes = vtkSmartPointer<vtkFillHolesFilter>::New();
        Holes -> SetInputData(Cell);
        Holes -> SetHoleSize(100.0);
        Holes -> Update();

        vtkSmartPointer<vtkSmoothPolyDataFilter> Smooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
        if (Holes->GetOutput()->GetNumberOfCells() > 0) {
            Smooth -> SetInputData(Holes->GetOutput());
        } else {
            Smooth -> SetInputData(Cell);
        }
        Smooth -> SetNumberOfIterations(20);
        Smooth -> SetRelaxationFactor(0.1);
        Smooth -> Update();

        vtkSmartPointer<vtkPolyDataNormals> Normals = vtkSmartPointer<vtkPolyDataNormals>::New();
        Normals -> SetInputData(Smooth->GetOutput());
        Normals -> ComputeCellNormalsOn();
        Normals -> ComputePointNormalsOn();
        Normals -> AutoOrientNormalsOn();
        Normals -> Update();

        Cell -> DeepCopy(Normals->GetOutput());
        Cell -> Modified();
    }

    void _Supernova::ScalePolyData(double _dxy, double _dz) {
        double r[3];
        vtkSmartPointer<vtkPoints> P = Rays -> GetPoints();
        for (vtkIdType i = 0; i < P->GetNumberOfPoints(); i++) {
            P -> GetPoint(i,r);
            r[0] *= _dxy; r[1] *= _dxy; r[2] *= _dz;
            P -> SetPoint(i,r);
        }
        P = Peaks -> GetPoints();
        for (vtkIdType i = 0; i < P->GetNumberOfPoints(); i++) {
            P -> GetPoint(i,r);
            r[0] *= _dxy; r[1] *= _dxy; r[2] *= _dz;
            P -> SetPoint(i,r);
        }
        P = Cell -> GetPoints();
        for (vtkIdType i = 0; i < P->GetNumberOfPoints(); i++) {
            P -> GetPoint(i,r);
            r[0] *= _dxy; r[1] *= _dxy; r[2] *= _dz;
            P -> SetPoint(i,r);
        }
    }
