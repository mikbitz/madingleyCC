#ifndef HERBIVORYTRACKER_H
#define HERBIVORYTRACKER_H
/** \file HerbivoryTracker.h
 * \brief the HerbivoryTracker header file
 */





//


//
//namespace Madingley
//{
/** \brief
//    /// Tracks the herbivory ecological process
//    /// </summary>
//     class HerbivoryTracker
//    {
/** \brief
The flow of mass between prey and predator
*/
//        private double[] _MassFlows;
/** \brief
Get and set the flow of mass between prey and predator
*/
//         double[] MassFlows
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
/** \brief
Get and set the missing data value to be used in the mass flows output
*/
//         double MissingValue { get { return _MissingValue; } set { _MissingValue = value; } }
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
The time steps to be run in the current simulation
*/
//        private float[] TimeSteps;
//
/** \brief
Set up the herbivory tracker
*/
@param numTimeSteps The total number of timesteps for this simulation 
@param cellIndices List of indices of active cells in the model grid 
@param massFlowsFilename Filename for outputs of the flows of mass between predators and prey 
@param cohortDefinitions The functional group definitions for cohorts in the model 
@param missingValue The missing value to be used in the output file 
@param outputFileSuffix The suffix to be applied to the output file 
@param outputPath The path to write the output file to 
@param trackerMassBins The mass bin handler containing the mass bins to be used for predation tracking 
//         HerbivoryTracker(unsigned numTimeSteps,
//            List<unsigned[]> cellIndices,
//            string massFlowsFilename,
//            FunctionalGroupDefinitions cohortDefinitions,
//            double missingValue,
//            string outputFileSuffix,
//            string outputPath, MassBinsHandler trackerMassBins)
//        {
//            // Assign the missing value
//            _MissingValue = missingValue;
//
//            // Get the mass bins to use for the predation tracker and the number of mass bins that this correpsonds to
//            _MassBins = trackerMassBins.GetSpecifiedMassBins();
//            _NumMassBins = trackerMassBins.NumMassBins;
//
//            // Initialise the array to hold data on mass flows between mass bins
//            _MassFlows = new double[_NumMassBins];
//
//            // Define the model time steps to be used in the output file
//            TimeSteps = new float[numTimeSteps];
//            for (int i = 1; i <= numTimeSteps; i++)
//            {
//                TimeSteps[i - 1] = i;
//            }
//
//            // Initialise the data converter
//            DataConverter = new ArraySDSConvert();
//
//            // Initialise the SDS object creator
//            SDSCreator = new CreateSDSObject();
//
//            // Create an SDS object to hold the predation tracker data
//            MassFlowsDataSet = SDSCreator.CreateSDS("netCDF", massFlowsFilename + outputFileSuffix, outputPath);
//
//            // Define the dimensions to be used in the predation tracker output file
//            string[] dimensions = { "Time step", "Herbivore mass bin" };
//
//            // Add the mass flow variable to the predation tracker
//            DataConverter.AddVariable(MassFlowsDataSet, "Log mass (g)", 2, dimensions, _MissingValue, TimeSteps, _MassBins);
//
//            
//        }
//
/** \brief
Record mass flow in an eating event
*/
@param timestep The current model time step 
@param herbivoreBiomass The individual body mass of the herbivore 
@param massFlow The amount of mass consumed in the predation event 
//         void RecordFlow(unsigned timestep, double herbivoreBiomass, double massFlow)
//        {
//
//            // Find the appropriate mass bin for the cohort
//            int HerbivoreMassBin = 0;
//            do
//            {
//                HerbivoreMassBin++;
//            } while (HerbivoreMassBin < (_MassBins.Length - 1) && herbivoreBiomass > _MassBins[HerbivoreMassBin]);
//
//            _MassFlows[HerbivoreMassBin] += massFlow;
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
//            string[] dimensions = { "Time step", "Herbivore mass bin" };
//
//            // Log all values of the mass flow
//            for (int i = 0; i < _NumMassBins; i++)
//            {
//                if (_MassFlows[i] > 0) _MassFlows[i] = Math.Log(_MassFlows[i]);
//                else _MassFlows[i] = _MissingValue;
//            }
//            // Add the mass flows data to the output file
//            DataConverter.VectorToSDS2D(_MassFlows, "Log mass (g)", dimensions, TimeSteps, _MassBins, _MissingValue, MassFlowsDataSet, timeStep);
//
//        }
//
/** \brief
Resets the mass flows data array
*/
//         void ResetHerbivoryTracker()
//        {
//            _MassFlows = new double[_NumMassBins];
//        }
//
/** \brief
Close the herbivory tracker
*/
//         void CloseStreams()
//        {
//            MassFlowsDataSet.Dispose();
//        }
//    }
//}
#endif
