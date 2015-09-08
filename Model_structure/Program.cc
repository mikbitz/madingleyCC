/// @todo complete the program file - this is just a case of setting the output directory
/// @todo complete the RunSimulations class - do we really need anything more? - scenarios should defined externally really...
/// @todo Complete the ScenarioParameterInitialisation class - again- are scenarios really internal to the code?
/// @todo Complete the MadingleyModelInitialisation class
/// @todo Complete the stopwatch class
/// @todo Complete the MadingleyModel class
/// @todo Complete the FunctionalGroupDefinitions class
/// @todo Complete the CohortMerge class
/// @todo Complete the ModelGrid class (a big one!)
/// @todo Complete the EcologyCrossGrid class
/// @todo Complete the Envirodata class
/// @todo Complete the MassBinsHanlder class
/// @todo Complete all teh dispersal implementations

/// @todo Test the UtilityFunctions class
/// @todo Test the Activity class
#include <iostream>
#include <math.h>

//Magindgley model entry point.
//Changes:-
//MB 6/8/2015
//Commented out all lines of code
//Created main function
//Created runsimulations header file
//Created Utilityfunction header file
//Created ScenarioParameterInitialisation header file
//Created MadingelyModelInitialzation header file
//Created stopwatch header file
//Created the MadingleyModel header file
//7/8/15
//Created FunctionalGroupDefinitions header file
//Created CohortMerge header file
//Created ModelGrid class beader
//Created the EcologyCrossGrid header
//Created the Envirodata header
//Created the MassBinHandler header
/**
\file Program.cc
\brief This is the main entry point for the madingley code
*/
/**
\brief Main program

Initialise and run the model

*/
#include <string>
#include <RunSimulations.h>
#include<vector>
#include <chrono>
#include <ctime>
#include <EcologyStock.h>
#include <IDispersalImplementation.h>

//
//
//namespace Madingley
//{

/** \brief Starts a model run or set of model runs
*/
///@todo Complete set up of output directory
        int main()
        {   

           // Write out model details to the console
            cout<<("Madingley model C++ v. 0.\n")<<endl;

           // Declare an instance of RunSimulations
            RunSimulations MakeSimulations = RunSimulations();
            //Outpu directory
            std::string OutputDir=".";
            
            std::time_t t = system_clock::to_time_t(high_resolution_clock::now());
            cout<<"Model Run started at "<<std::ctime(&t)<<endl;



           // Declare an instance of ScenarioParameterInitialisation to read in the parameters for this model run or set of runs
            ScenarioParameterInitialisation Scenarios = ScenarioParameterInitialisation("Scenarios.csv", OutputDir);

           // Run the desired simulation or batch of simulations
            MakeSimulations.RunAllSimulations("EcosystemModelInitialisation.csv", Scenarios, OutputDir);
       
            return 0;
    }
//}