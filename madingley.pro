 HEADERS       = ./Model_structure/MadingleyModelInitialisation.h \
                 ./Model_structure/MadingleyModel.h \
                 ./Model_structure/ModelGrid.h \
                 ./Utility_classes/Properties.h \
                 ./Utility_classes/Stopwatch.h 
 SOURCES       = ./Model_structure/Program.cc \
                 ./Model_structure/Stock.cc \
                 ./Model_structure/Cohort.cc \
                 ../NetCDFIO/src/Convertor.cpp \
                 ../NetCDFIO/src/DataGrid.cpp \
                 ../NetCDFIO/src/FileReader.cpp \
                 ../NetCDFIO/src/NcGridCell.cpp \
                 ../NetCDFIO/src/Logger.cpp 
 RESOURCES     = 
 INCLUDEPATH += ../Madingley/Model_structure ../Madingley/Utility_classes ../Madingley/Ecological_processes_cohorts ../Madingley/Ecological_processes_cohorts/Technical_code ../Madingley/Output_and_tracking 
 INCLUDEPATH += ../Madingley/Ecological_processes_cohorts/Dispersal_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Eating_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Metabolism_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Mortality_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Reproduction_implementations/Technical_code
 INCLUDEPATH += ../NetCDFIO/src /usr/local/netcdfcxx-4.2.1/include /usr/local/netcdfc-4.3.3.1/include
 DEPENDPATH += ../Madingley/Model_structure ../Madingley/Utility_classes ../Madingley/Ecological_processes_cohorts ../Madingley/Ecological_processes_cohorts/Technical_code ../Madingley/Output_and_tracking
 DEPENDPATH += ../Madingley/Ecological_processes_cohorts/Dispersal_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Eating_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Metabolism_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Mortality_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Reproduction_implementations/Technical_code
 
 QMAKE_CXXFLAGS += -std=c++11 -w
 QMAKE_LIBDIR += /usr/local/netcdfcxx-4.2.1/lib64 /usr/local/netcdfc-4.3.3.1/lib64 /usr/local/hdf-1.8.15/lib64 /usr/local/zlib-1.2.8/lib64
 LIBS += -lnetcdf_c++4 -lnetcdf -lhdf5_hl -lhdf5 -lcurl -lz -ldl
 # install
 target.path = madingley 
 sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS  images
 sources.path = .:./Model_structure
 INSTALLS += target sources


