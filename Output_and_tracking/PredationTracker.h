#ifndef PREDATIONTRACKER_H
#define PREDATIONTRACKER_H
/** \file PredationTracker.h
 * \brief the PredationTracker header file
 */





//


//
//namespace Madingley
//{
/** \brief
//    /// Tracks the predation ecological process
//    /// </summary>
//     class PredationTracker
//    {
/** \brief
The flow of mass between prey and predator
*/
//        private double[,] _MassFlows;
/** \brief
Get and set the flow of mass between prey and predator
*/
//         double[,] MassFlows
//        {
//            get { return _MassFlows; }
//            set { _MassFlows = value; }
//        }
//
//
/** \brief
Vector of mass bins to be used in the predation tracker
*/
//        private float[] _MassBins;
//
/** \brief
The number of mass bins to track predation for
*/
//        private int _NumMassBins;
//
//
/** \brief
Missing data value to be used in the mass flows output
*/
//        private double _MissingValue;
//    	/// <summary>
Get and set the missing data value to be used in the mass flows output
//    	/// </summary>
//         double MissingValue {get { return _MissingValue;} set { _MissingValue = value;} }
//	
//
/** \brief
Dataset to output the Massflows data
*/
//        private DataSet MassFlowsDataSet;
//
//
/** \brief
An instance of the class to convert data between arrays and SDS objects
*/
//        private ArraySDSConvert DataConverter;
//
/** \brief
Instance of the class to create SDS objects
*/
//        private CreateSDSObject SDSCreator;
//
/** \brief
Set up the predation tracker
*/_
@param numTimeSteps The total number of timesteps for this simulation 
@param cellIndices List of indices of active cells in the model grid 
@param massFlowsFilename Filename for outputs of the flows of mass between predators and prey 
@param cohortDefinitions The functional group definitions for cohorts in the model 
@param missingValue The missing value to be used in the output file 
@param outputFileSuffix The suffix to be applied to the output file 
@param outputPath The path to write the output file to 
@param trackerMassBins The mass bin handler containing the mass bins to be used for predation tracking 
//         PredationTracker(unsigned numTimeSteps,
//            List<unsigned[]> cellIndices, 
//            string massFlowsFilename, 
//            FunctionalGroupDefinitions cohortDefinitions, 
//            double missingValue,
//            string outputFileSuffix,
//            string outputPath, MassBinsHandler trackerMassBins, int cellIndex)
//        {
//            // Assign the missing value
//            _MissingValue = missingValue;
//
//            // Get the mass bins to use for the predation tracker and the number of mass bins that this correpsonds to
//            _MassBins = trackerMassBins.GetSpecifiedMassBins();
//            _NumMassBins = trackerMassBins.NumMassBins;
//
//            // Initialise the array to hold data on mass flows between mass bins
//            _MassFlows = new double[_NumMassBins, _NumMassBins];
//            
//            // Define the model time steps to be used in the output file
//            float[] TimeSteps = new float[numTimeSteps];
//            for (int i = 1; i <= numTimeSteps; i++)
//            {
//                TimeSteps[i-1] = i;
//            }
//
//            // Initialise the data converter
//            DataConverter = new ArraySDSConvert();
//
//            // Initialise the SDS object creator
//            SDSCreator = new CreateSDSObject();
//
//            // Create an SDS object to hold the predation tracker data
//            MassFlowsDataSet = SDSCreator.CreateSDS("netCDF", massFlowsFilename + outputFileSuffix + "_Cell" + cellIndex, outputPath);
//
//            // Define the dimensions to be used in the predation tracker output file
//            string[] dimensions = { "Predator mass bin", "Prey mass bin", "Time steps" };
//
//            // Add the mass flow variable to the predation tracker
//            DataConverter.AddVariable(MassFlowsDataSet, "Log mass (g)", 3, dimensions, _MissingValue, _MassBins, _MassBins, TimeSteps);    
//        }
//
/** \brief
Record mass flow in an eating event
*/
@param timestep The current model time step 
@param preyBiomass The individual body mass of the prey 
@param predatorBiomass The individual body mass of the predator 
@param massFlow The amount of mass consumed in the predation event 
//         void RecordFlow(unsigned timestep, double preyBiomass, double predatorBiomass, double massFlow)
//        {
//            
//            // Find the appropriate mass bin for the cohort
//            int PredatorMassBin = 0;
//            do
//            {
//                PredatorMassBin++;
//            } while (PredatorMassBin < (_MassBins.Length - 1) && predatorBiomass > _MassBins[PredatorMassBin]);
//
//            // Find the appropriate mass bin for the cohort
//            int PreyMassBin = 0;
//            do
//            {
//                PreyMassBin++;
//            } while (PreyMassBin < (_MassBins.Length - 1) && preyBiomass > _MassBins[PreyMassBin]);
//
//            _MassFlows[PredatorMassBin,PreyMassBin] += massFlow;
//
//        }
//
/** \brief
Add the mass flows from the current timestep to the dataset
*/
@param timeStep the current timestep 
//         void AddTimestepFlows(int timeStep)
//        {
//            // Define the dimensions of the output data
//            string[] dimensions = { "Predator mass bin", "Prey mass bin", "Time steps" };
//
//            // Log all values of the mass flow
//            for (int i = 0; i < _NumMassBins; i++)
//            {
//                for (int j = 0; j < _NumMassBins; j++)
//                {
//                    if (_MassFlows[i, j] > 0) _MassFlows[i, j] = Math.Log(_MassFlows[i, j]);
//                    else _MassFlows[i, j] = _MissingValue;
//                }
//            }
//            // Add the mass flows data to the output file
//            DataConverter.Array2DToSDS3D(_MassFlows, "Log mass (g)", dimensions, timeStep, _MissingValue, MassFlowsDataSet);
//
//        }
//
/** \brief
Resets the mass flows data array
*/
//         void ResetPredationTracker()
//        {
//            _MassFlows = new double[_NumMassBins, _NumMassBins];
//        }
//
//
/** \brief
Close the predation tracker
*/
//         void CloseStreams()
//        {
//            MassFlowsDataSet.Dispose();
//        }
//    }
//}
#endif
