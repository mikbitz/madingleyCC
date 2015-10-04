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
    class RunSimulations
    {
    public:
/**       
\brief
Runs the specified number of simulations for each of the specified scenarios

@param initialisationFilename Filename of the file from which to read initialisation information
@param scenarios Contains scenario information for this set of simulations
@param outputPath The path to which outputs should be written
*/
         void RunAllSimulations(string initialisationFilename, ScenarioParameterInitialisation scenarios, string outputPath)
        {
//            // Declare an instance of the class for initializing the Madingley model
            MadingleyModelInitialisation InitialiseMadingley = MadingleyModelInitialisation(initialisationFilename, outputPath);
//            // Specify the output path in this instance
            InitialiseMadingley.OutputPath = outputPath;
//
//           
//
//            // List to hold the names of the scenarios to run
            vector<string> ScenarioNames;
//            // String variable to hold the index suffix to apply to output files for a given simulation
            string OutputFilesSuffix;
//
//            // Loop over scenario names and add the name of the scenario to the list of scenarion names
//            foreach (string scenarioName in scenarios.scenarioParameters.Keys)
//            {
//                ScenarioNames.Add(scenarioName);
//            }
//            
//            // Check whether there is only one simulation to run
//            if (scenarios.scenarioNumber == 1 && scenarios.scenarioSimulationsNumber[scenarios.scenarioNumber-1] == 1)
//            {
//                // For a single simulation
//
//                // Set-up the suffix for the output files
                OutputFilesSuffix = "_";
//                
//                // Loop over the parameters for this scenario
//                for (int i = 0; i < ScenarioNames.Count; i++)
//                {
//                    // Add the parameter information to the suffix for this simulation
//                    OutputFilesSuffix += scenarios.scenarioParameters[ScenarioNames[i]].ElementAt(0) + "_";
//                }
//                // Add a zero index to the end of the suffix
//                OutputFilesSuffix += "0";
//
//                //Run the simulation
                RunSimulation(scenarios, 0, InitialiseMadingley, OutputFilesSuffix);
//                
//            }
//            else
//            {
//
//                if (InitialiseMadingley.RunSimulationsInParallel && InitialiseMadingley.InitialisationFileStrings.ContainsKey("Locations"))
//                {
//                    // Loop over specified scenarios iteratively
//                    for (int ScenarioIndex = 0; ScenarioIndex < scenarios.scenarioNumber; ScenarioIndex++)
//                    {
//                        //Create an array of new MadingleyModel instances for simulations under this scenario combination
//                        MadingleyModel[] MadingleyEcosystemModels = new MadingleyModel[scenarios.scenarioSimulationsNumber[ScenarioIndex]];
//                        
//                        for (int simulation = 0; simulation < scenarios.scenarioSimulationsNumber[ScenarioIndex]; simulation++)
//                        {
//                            // Set up the suffix for the output files
//                            OutputFilesSuffix = "_";
//
//                            // Loop over parameters for this scenario
//                            for (int i = 0; i < ScenarioNames.Count; i++)
//                            {
//                                // Add the parameter information to the suffix for this simulation
//                                OutputFilesSuffix += scenarios.scenarioParameters[ScenarioNames[i]].ElementAt(ScenarioIndex) + "_";
//                            }
//                            // Add the simulation index number to the suffix
//                            OutputFilesSuffix += simulation.ToString();
//
//                            // Initialize the instance of MadingleyModel
//                            MadingleyEcosystemModels[simulation] = new MadingleyModel(InitialiseMadingley, scenarios, ScenarioIndex, OutputFilesSuffix,
//                                InitialiseMadingley.GlobalModelTimeStepUnit);
//                        }
//
//                        // Loop over the specified number of simulations for each scenario
//                        //for (int simulation = 0; simulation<  scenarios.scenarioSimulationsNumber[ScenarioIndex]; simulation++)
//                        Parallel.For(0, scenarios.scenarioSimulationsNumber[ScenarioIndex], simulation =>
//                        {
//                            // Declare and start a timer
//                            StopWatch s = new StopWatch();
//                            s.Start();
//
//                            // Run the simulation
//                            MadingleyEcosystemModels[simulation].RunMadingley(InitialiseMadingley);
//
//                            // Stop the timer and write out the time taken to run this simulation
//                            s.Stop();
//                            Console.WriteLine("Model run finished");
//                            Console.WriteLine("Total elapsed time was {0} seconds", s.GetElapsedTimeSecs());
//
//                        });
//                    }
//                }
//                else
//                {
//                    //Run simulations sequentially
//
//                    // Loop over specified scenarios
//                    for (int ScenarioIndex = 0; ScenarioIndex < scenarios.scenarioNumber; ScenarioIndex++)
//                    {
//                        // Loop over the specified number of simulations for each scenario
//                        for (int simulation = 0; simulation < scenarios.scenarioSimulationsNumber[ScenarioIndex]; simulation++)
//                        {
//                            // Set up the suffix for the output files
//                            OutputFilesSuffix = "_";
//
//                            // Loop over parameters for this scenario
//                            for (int i = 0; i < ScenarioNames.Count; i++)
//                            {
//                                // Add the parameter information to the suffix for this simulation
//                                OutputFilesSuffix += scenarios.scenarioParameters[ScenarioNames[i]].ElementAt(ScenarioIndex) + "_";
//                            }
//                            // Add the simulation index number to the suffix
//                            OutputFilesSuffix += simulation.ToString();
//
//                            // Run the current simulation
//                            RunSimulation(scenarios, ScenarioIndex, InitialiseMadingley, OutputFilesSuffix);
//                        }
//                    }
//                }
//               
//            }
//
        }
//
/** \brief Runs a single simulation of the Madingley model 
@param scenarios Parameter information and simulation number for all scenarios to be run
@param scenarioIndex The index of the scenario to be run in this simulation
@param initialiseMadingley Model initialization information for all simulations
@param outputFileSuffix Suffix to be applied to the names of files written out by this simulation */
void RunSimulation(ScenarioParameterInitialisation scenarios, int scenarioIndex, MadingleyModelInitialisation initialiseMadingley, 
            string outputFileSuffix)
        {
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
            cout<<"Model run finished"<<endl;
            cout<<"Total elapsed time was "<< s.GetElapsedTimeSecs() << " seconds "<<endl;
        }

};
#endif