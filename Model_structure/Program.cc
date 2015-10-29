#include <iostream>
#include <math.h>
#include <fenv.h>
#include <cstdlib>


//Magindgley model entry point.
//Changes from original PLOS Biology version:-
//MB 6/8/2015 - 4/10/2015 - converted C# to C++

/** \file Program.cc \brief This is the main entry point for the madingley code */
/**
\brief Main program

Initialise and run the model

 */
#include <MadingleyModel.h>
//----------------------------------------------------------------------------------------------

/** \brief Starts a model run or set of model runs */
///@todo Complete set up of output directory

int main() {
    //this line enables the gdb debugger to catch Nan or floating point problems
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    // Write out model details to the console
    cout << ("Madingley model C++ v. 0.\n") << endl;

    //Output directory
    std::string OutputPath = ".";
    //model run parameters
    std::string initialisationFileName = "input/Model setup/EcosystemModelInitialisation.csv";

    std::time_t t = system_clock::to_time_t(high_resolution_clock::now());
    cout << "Model Run started at " << std::ctime(&t) << endl;

    // Initialise the model
    // Declare an instance of the class that runs a Madingley model simulation
    MadingleyModel MadingleyEcosystemModel(initialisationFileName,OutputPath);

    // Declare and start a timer
    StopWatch s;
    s.Start();

    // Run the simulation
    MadingleyEcosystemModel.RunMadingley();

    // Stop the timer and write out the time taken to run this simulation
    s.Stop();
    cout << "Model run finished" << endl;
    cout << "Total elapsed time was " << s.GetElapsedTimeSecs() << " seconds " << endl;
    return 0;
}
//----------------------------------------------------------------------------------------------
