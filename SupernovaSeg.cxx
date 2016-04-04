
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

    #ifdef DEBUG
        printf("SuperNova Segmentation V1.0 [DEBUG mode]\n");
        printf("File name: %s\n",DataBase->GetFullCellName().c_str());
    #endif

    vtkSmartPointer<vtkTIFFReader> TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
    TIFFReader -> SetOrientationType(TIFF_ORIENTATION_READER);
    TIFFReader -> SetFileName(DataBase->GetFullCellName().c_str());
    TIFFReader -> Update();

    vtkSmartPointer<vtkImageData> ImageData = TIFFReader -> GetOutput();
    vtkSmartPointer<vtkImageData> ClipImage = vtkSmartPointer<vtkImageData>::New();

    int *Dim = ImageData -> GetDimensions();

    #ifdef DEBUG
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

        Supernova -> SetScaleFactor(DataBase->GetDxy(),DataBase->GetDz());
        
        Supernova -> Initialize();
        
        Supernova -> Probe(ImageData);
        
        if (DataBase->GetR1(id)>20) {

            Supernova -> ApplyLimits(DataBase->GetR1(id),false);

        } else {

            #ifdef DEBUG
                printf("Small cell detected. Assuming spherical shape...\n");
            #endif

            Supernova -> ApplyLimits(DataBase->GetR1(id),true);

        }
        
        Supernova -> Segmentation();

        Supernova -> ClipImageData(DataBase->GetFullMitoName().c_str(),ImageData,ClipImage,DataBase->GetId(id));
        
        SaveImageData(DataBase->MakeGenericFileName(DataBase->GetId(id),"-mitovolume",".vtk").c_str(),ClipImage);

        Supernova -> ScalePolyData(DataBase->GetDxy(),DataBase->GetDz());

        Supernova -> Save(DataBase->MakeGenericFileName(DataBase->GetId(id),"-cellsurface",".vtk").c_str(),DataBase->MakeGenericFileName(DataBase->GetId(id),"-celloutersurface",".vtk").c_str());

        Supernova -> SaveMassProperties(DataBase->MakeGenericFileName(DataBase->GetId(id),"",".supernovaseg").c_str());

    }   

    return EXIT_SUCCESS;
}

int ScanFolderForThisExtension(const char _root[], const char ext[], std::vector<std::string> *List) {
    DIR *dir;
    int ext_p;
    struct dirent *ent;
    std::string _dir_name;
    if ((dir = opendir (_root)) != NULL) {
      while ((ent = readdir (dir)) != NULL) {
        _dir_name = std::string(ent->d_name);
        ext_p = (int)_dir_name.find(std::string(ext));
        if (ext_p > 0) {
            #ifdef DEBUG
                printf("File found: %s\n",_dir_name.c_str());
            #endif
            List -> push_back(std::string(_root)+_dir_name.substr(0,ext_p));
        }
      }
      closedir (dir);
    } else {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

/* ================================================================
   MAIN ROUTINE
=================================================================*/

int main(int argc, char *argv[]) {     
    
    std::string Path;
    _database *DataBase = new _database;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            Path = std::string(argv[i+1]);
            if (Path.back() != '/')
                Path = Path + '/';
        }
        if (!strcmp(argv[i],"-check")) {
            DataBase->SetCheckModeOn();
        }
    }

    std::vector<std::string> Files;
    ScanFolderForThisExtension(Path.c_str(),".centers",&Files);

    for (int i = 0; i < Files.size(); i++) {
        DataBase -> PopulateFromFile(Files[i]+".centers");
        RunSuperNova(DataBase);
        if ( DataBase->CheckMode() ) {
            #ifdef DEBUG
                printf("Running check mode...\n");
            #endif
            break;
        }
    }

    return EXIT_SUCCESS;
}
