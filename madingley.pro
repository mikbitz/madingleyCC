 HEADERS       = ./Model_structure/RunSimulations.h \
                 ./Model_structure/ScenarioParameterInitialisation.h \
                 ./Model_structure/MadingleyModelInitialisation.h \
                 ./Utility_classes/Properties.h \
                 ./Utility_classes/Stopwatch.h
 SOURCES       = ./Model_structure/Program.cc 
 RESOURCES     = 
 INCLUDEPATH += ../Madingley/Model_structure ../Madingley/Utility_classes ../Madingley/Ecological_processes_cohorts ../Madingley/Ecological_processes_cohorts/Technical_code ../Madingley/Output_and_tracking 
 INCLUDEPATH += ../Madingley/Ecological_processes_cohorts/Dispersal_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Eating_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Metabolism_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Mortality_implementations/Technical_code ../Madingley/Ecological_processes_cohorts/Reproduction_implementations/Technical_code

 # install
 target.path = madingley 
 sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS  images
 sources.path = .:./Model_structure
 INSTALLS += target sources


