#ifndef RUNSIMULATIONS_H
#define RUNSIMULATIONS_H
#include <string>
#include <ScenarioParameterInitialisation.h>
#include <MadingleyModelInitialisation.h>
#include <MadingleyModel.h>
#include <Stopwatch.h>

/**
\file RunSimulations.h
\brief The Runs simulations header file
 */

using namespace std;

/**
\class RunSimulations
\brief Runs simulations of the Madingley model

 */
class RunSimulations {
public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    /** \brief Runs the specified number of simulations for each of the specified scenarios

    @param initialisationFilename Filename of the file from which to read initialisation information
    @param scenarios Contains scenario information for this set of simulations
    @param outputPath The path to which outputs should be written
     */
    void RunAllSimulations(string initialisationFilename, ScenarioParameterInitialisation scenarios, string outputPath) {
        // Declare an instance of the class for initializing the Madingley model
        MadingleyModelInitialisation InitialiseMadingley = MadingleyModelInitialisation(initialisationFilename, outputPath);
        // Specify the output path in this instance
        InitialiseMadingley.OutputPath = outputPath;

        // List to hold the names of the scenarios to run
        vector<string> ScenarioNames;
        // String variable to hold the index suffix to apply to output files for a given simulation
        string OutputFilesSuffix;

        // Set-up the suffix for the output files
        OutputFilesSuffix = "_";

        //Run the simulation
        RunSimulation(scenarios, 0, InitialiseMadingley, OutputFilesSuffix);
    }
    
    //----------------------------------------------------------------------------------------------
    /** \brief Runs a single simulation of the Madingley model 
    @param scenarios Parameter information and simulation number for all scenarios to be run
    @param scenarioIndex The index of the scenario to be run in this simulation
    @param initialiseMadingley Model initialization information for all simulations
    @param outputFileSuffix Suffix to be applied to the names of files written out by this simulation */
    void RunSimulation(ScenarioParameterInitialisation scenarios, int scenarioIndex, MadingleyModelInitialisation initialiseMadingley,
            string outputFileSuffix) {
        // Declare an instance of the class that runs a Madingley model simulation
        MadingleyModel MadingleyEcosystemModel(initialiseMadingley, scenarios, scenarioIndex, outputFileSuffix,
                initialiseMadingley.GlobalModelTimeStepUnit);

        // Declare and start a timer
        StopWatch s;
        s.Start();

        // Run the simulation
        MadingleyEcosystemModel.RunMadingley(initialiseMadingley);

        // Stop the timer and write out the time taken to run this simulation
        s.Stop();
        cout << "Model run finished" << endl;
        cout << "Total elapsed time was " << s.GetElapsedTimeSecs() << " seconds " << endl;
    }
    //----------------------------------------------------------------------------------------------

};
#endif