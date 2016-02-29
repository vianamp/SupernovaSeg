
#include "_database.h"
#include "_Supernova.h"

int SaveImageData(const char FileName[], vtkImageData *Image) {
    vtkSmartPointer<vtkStructuredPointsWriter> W = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    W -> SetInputData(Image);
    W -> SetFileName(FileName);
    W -> Write();
    return EXIT_SUCCESS;
}

int RunSuperNova(_database *DataBase) {

    vtkSmartPointer<vtkTIFFReader> TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
    TIFFReader -> SetOrientationType(1);
    TIFFReader -> SetFileName(DataBase->GetFullSurfaceName().c_str());
    TIFFReader -> Update();

    vtkSmartPointer<vtkImageData> ImageData = TIFFReader -> GetOutput();
    vtkSmartPointer<vtkImageData> ClipImage = vtkSmartPointer<vtkImageData>::New();

    int *Dim = ImageData -> GetDimensions();

    #ifdef DEBUG
        printf("SuperNova Segmentation V1.0 [DEBUG mode]\n");
        printf("File name: %s\n",DataBase->GetFullSurfaceName().c_str());
        printf("Volume dimensions: %dx%dx%d\n",Dim[0],Dim[1],Dim[2]);
    #endif

    int xo, yo, zo, tp;
    for (int id = 0; id < DataBase->GetNumberOfCenters(); id++) {

        xo = DataBase -> GetX(id);
        yo = DataBase -> GetY(id);
        zo = DataBase -> GetZ(id);
        tp = DataBase -> GetType(id);

        printf("Running Center [%d (%d)]: %d, %d, %d\n",id,tp,xo,yo,zo);

        _Supernova *Supernova = new _Supernova();
        
        Supernova -> SetCenter(xo,yo,zo);
        
        Supernova -> Initialize();
        
        Supernova -> Probe(ImageData);
        
        if ( tp == 2 ) {
            Supernova -> ApplyLimits(DataBase->GetR1(id),false);
        } else if ( tp == 3 ) {
            Supernova -> ApplyLimits(DataBase->GetR1(id),true);
        }
        
        Supernova -> Segmentation(id);

        Supernova -> ClipImageData(DataBase->GetFullMitoName().c_str(),ImageData,ClipImage);

        Supernova -> AjustCoordinates(Dim[1],DataBase->GetDxy(),DataBase->GetDz());
        
        Supernova -> Save(DataBase->MakeVTKFileName(id,"-surface").c_str(),id);

        SaveImageData(DataBase->MakeVTKFileName(id,"-volume").c_str(),ClipImage);

    }

    return EXIT_SUCCESS;
}

/* ================================================================
   MAIN ROUTINE
=================================================================*/

int main(int argc, char *argv[]) {     
    
    _database *DataBase = new _database;

    int i;
    bool _pts = 0;
    char _prefix[64];
    char _impath[128];
    sprintf(_impath,"");

    for (i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            sprintf(_impath,"%s",argv[i+1]);
        }
    }

    char _cmd[256];
    sprintf(_cmd,"ls %s*.centers > %ssupernovaseg.files",_impath,_impath);
    system(_cmd);

    char _tifffilename[256];
    char _tifflistpath[128];
    sprintf(_tifflistpath,"%ssupernovaseg.files",_impath);
    FILE *f = fopen(_tifflistpath,"r");

    while (fgets(_tifffilename,256, f) != NULL) {
        _tifffilename[strcspn(_tifffilename, "\n" )] = '\0';
        DataBase -> PopulateFromFile(_tifffilename);
        RunSuperNova(DataBase);
        return 0;
    }

    fclose(f);

    return EXIT_SUCCESS;
}
