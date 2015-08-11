#ifndef SCENARIOPARAMETERINITIALISATION_H
#define SCENARIOPARAMETERINITIALISATION_H
#include <string>
#include <Properties.h>
using namespace std;
/**
\file ScenarioParameterInitialisation.h
\brief The ScenarioParameterInitialisation header file
*/



//
//using Microsoft.Research.Science.Data;
//
//namespace Madingley
//{
/**
\class  ScenarioParameterInitialisation
\brief Reads the file specifying which scenarios will be run, and stores this information
*/
class ScenarioParameterInitialisation{
public:
    /** The number of scenarios to be run*/
    IntProperty scenarioNumber; 
//
//        /// <summary>
//        /// String parameters for the different scenarios to run
//        /// </summary>
//        private SortedList<string, string[]> _scenarioParameters;
//        /// <summary>
//        /// Get the string parameters for the different scenarios to run
//        /// </summary>
//        public SortedList<string, string[]> scenarioParameters
//        { get { return _scenarioParameters; } }
//
//        /// <summary>
//        /// The number of simulations to run for each of the scenarios
//        /// </summary>
//        private List<int> _scenarioSimulationsNumber;
//        /// <summary>
//        /// Get the number of simulations to run for each of the scenarios
//        /// </summary>
//        public List<int> scenarioSimulationsNumber
//        {
//            get { return _scenarioSimulationsNumber; }
//            set { _scenarioSimulationsNumber = value; }
//        }
public:        
/** \brief
Constructor for ScenarioParameterInitialisation: reads in scenario parameters from a specified file

<param name="scenarioParameterFile">The name of the scenario parameters file, which must be in the 'Model setup' directory</param>
<param name="outputPath">The directory to write output files to</param>
*/
ScenarioParameterInitialisation(string scenarioParameterFile, string outputPath)
        {
//            Console.WriteLine("Reading scenario parameters file...\n");
//
//            // Construct file name
//            string FileString = "msds:csv?file=input/Model setup/" + scenarioParameterFile + "&openMode=readOnly";
//
//            //Copy the scenarioParameterFile to the output directory
//            System.IO.File.Copy("input/Model setup/" + scenarioParameterFile, outputPath + scenarioParameterFile, true);
//
//            // Read in the data
//            DataSet InternalData = DataSet.Open(FileString);
//
//            // Get the number of scenarios to be run based on the number of lines in the first variable in the input file
//            _scenarioNumber = InternalData.Variables[0].GetData().Length;
//
//            // Intialise sorted lists for parameter combinations and simulations numbers
//            _scenarioParameters = new SortedList<string, string[]>();
//            _scenarioSimulationsNumber = new List<int>();
//
//            // Temporary vector to hold parameter information
//            string[] TempExtractionParameters = new string[_scenarioNumber];
//
//            // Loop over parameters in the scenarios file
//            foreach (Variable v in InternalData.Variables)
//            {
//                //Get the name of the variable currently referenced in the dataset
//                string HeaderName = v.Name;
//                //Copy the values for this variable into an array
//                var TempValues = v.GetData();
//
//                // Switch based on the header name for the scenario parameter
//                switch (HeaderName.ToLower())
//                {
//                    case "human npp extraction":
//                        // Loop over scenarios and extract the parameter values for each scenario
//                        for (int i = 0; i < _scenarioNumber; i++)
//                        {
//                            TempExtractionParameters[i] = TempValues.GetValue(i).ToString();
//                        }
//                        break;
//                    case "number of simulations":
//                        // Loop over scenarios and extract the number of simulations to run for each scenario
//                        for (int i = 0; i < _scenarioNumber; i++)
//                        {
//                            _scenarioSimulationsNumber.Add(Convert.ToInt32(TempValues.GetValue(i)));
//                        }
//                        break;
//                }
//            }
//
//            // Add the extracted parameter values to the sorted list of parameters
//            _scenarioParameters.Add("Human NPP Extraction", TempExtractionParameters);
//
//
//
//        }
//
    }
};
#endif