#ifndef ECOLOGYCOHORT_H
#define ECOLOGYCOHORT_H
/** \file EcologyCohort.h
 * \brief the EcologyCohort header file
 */








//
//namespace Madingley
//{
/** \brief
//    /// A class to specify, initalise and run ecological processes pertaining to cohorts
//    /// </summary>
//     class EcologyCohort
//    {
//        # region Properties and Fields
//
/** \brief
A vector of stopwatch objects for timing the ecological processes
*/
//         StopWatch[] s2;
//               
/** \brief
A sorted list of formulations of metabolism
*/
//        private map<string, IEcologicalProcessWithinGridCell> _MetabolismFormulations;
/** \brief
Get the sorted list of metabolism formulations
*/
//         map<string, IEcologicalProcessWithinGridCell> MetabolismFormulations
//	    {
//		    get { return _MetabolismFormulations;}
//        }
//
/** \brief
A sorted list of formulations of eating
*/
//        private map<string, IEcologicalProcessWithinGridCell> _EatingFormulations;
/** \brief
Get the sorted list of eating formulations
*/
//         map<string, IEcologicalProcessWithinGridCell> EatingFormulations
//	    {
//		    get { return _EatingFormulations;}
//        }
//	
/** \brief
A sorted list of formulations of mortality
*/
//        private map<string, IEcologicalProcessWithinGridCell> _MortalityFormulations;
/** \brief
Get the sorted list of mortality formulations
*/
//         map<string, IEcologicalProcessWithinGridCell> MortalityFormulations
//	    {
//		    get { return _MortalityFormulations;}
//        }
//       
/** \brief
A sorted list of formulations of reproduction
*/
//        private map<string, IEcologicalProcessWithinGridCell> _ReproductionFormulations;
/** \brief
Get the sorted list of reproduction formulations
*/
//         map<string, IEcologicalProcessWithinGridCell> Reproductions
//	    {
//		    get { return _ReproductionFormulations;}
//        }
//
/** \brief
An instance of apply ecology
*/
//        ApplyEcology ApplyEcologicalProcessResults;
//
//
//        # endregion
//
/** \brief
Initalise the ecological processes
*/
//         void InitializeEcology(double cellArea, string globalModelTimeStepUnit, Boolean drawRandomly)
//        {
//            // Initialise eating formulations
//            _EatingFormulations = new map<string, IEcologicalProcessWithinGridCell>();
//            // Declare and attach eating formulations
//            Eating EatingFormulation = new Eating(cellArea, globalModelTimeStepUnit);
//            _EatingFormulations.Add("Basic eating", EatingFormulation);
//
//            // Initialise metabolism formulations
//            _MetabolismFormulations = new map<string, IEcologicalProcessWithinGridCell>();
//            // Declare and attach metabolism formulations
//            Metabolism MetabolismFormulation = new Metabolism(globalModelTimeStepUnit);
//            _MetabolismFormulations.Add("Basic metabolism", MetabolismFormulation);
//
//            // Initialise mortality formulations
//            _MortalityFormulations = new map<string, IEcologicalProcessWithinGridCell>();
//            // Declare and attach mortality formulations
//            Mortality MortalityFormulation = new Mortality(globalModelTimeStepUnit);
//            _MortalityFormulations.Add("Basic mortality", MortalityFormulation);
//
//            // Initialise reproduction formulations
//            _ReproductionFormulations = new map<string, IEcologicalProcessWithinGridCell>();
//            // Declare and attach mortality formulations
//            Reproduction ReproductionFormulation = new Reproduction(globalModelTimeStepUnit, drawRandomly);
//            _ReproductionFormulations.Add("Basic reproduction", ReproductionFormulation);
//
//            // Initialise apply ecology
//            ApplyEcologicalProcessResults = new ApplyEcology();
//
//
//        }
//
/** \brief
Run ecological processes that operate on cohorts within a single grid cell
*/
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param actingCohort The acting cohort 
@param cellEnvironment The environment in the current grid cell 
@param deltas A sorted list of deltas to track changes in abundances and biomasses during the ecological processes 
@param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
@param madingleyStockDefinitions The definitions for stock functional groups in the model 
@param currentTimestep The current model time step 
@param trackProcesses An instance of the process tracker 
@param partial Thread-locked local variables 
@param specificLocations Whether the model is being run for specific locations 
@param outputDetail The level of output detail being used for this model run 
@param currentMonth The current model month 
//         void RunWithinCellEcology(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, int[] actingCohort, 
//            map<string, double[]> cellEnvironment, Dictionary<string, Dictionary<string, double>> deltas, FunctionalGroupDefinitions 
//            madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions, unsigned currentTimestep, ProcessTracker trackProcesses, 
//            ref ThreadLockedParallelVariables partial, Boolean specificLocations,string outputDetail, unsigned currentMonth, MadingleyModelInitialisation initialisation)
//        {
//
//            // RUN EATING
//            _EatingFormulations["Basic eating"].RunEcologicalProcess(gridCellCohorts, gridCellStocks, actingCohort, cellEnvironment,
//                deltas, madingleyCohortDefinitions, madingleyStockDefinitions, currentTimestep, trackProcesses, ref partial,
//                specificLocations, outputDetail, currentMonth, initialisation);
//
//            
//            // RUN METABOLISM - THIS TIME TAKE THE METABOLIC LOSS TAKING INTO ACCOUNT WHAT HAS BEEN INGESTED THROUGH EATING
//            _MetabolismFormulations["Basic metabolism"].RunEcologicalProcess(gridCellCohorts, gridCellStocks, actingCohort,
//                cellEnvironment, deltas, madingleyCohortDefinitions, madingleyStockDefinitions, currentTimestep, trackProcesses, ref partial,
//                specificLocations, outputDetail, currentMonth, initialisation);
//              
//           
//            // RUN REPRODUCTION - TAKING INTO ACCOUNT NET BIOMASS CHANGES RESULTING FROM EATING AND METABOLISING
//            _ReproductionFormulations["Basic reproduction"].RunEcologicalProcess(gridCellCohorts, gridCellStocks, actingCohort,
//                cellEnvironment, deltas, madingleyCohortDefinitions, madingleyStockDefinitions, currentTimestep, trackProcesses, ref partial,
//                specificLocations, outputDetail, currentMonth, initialisation);
//            
//              
//            // RUN MORTALITY - TAKING INTO ACCOUNT NET BIOMASS CHANGES RESULTING FROM EATING, METABOLISM AND REPRODUCTION
//            _MortalityFormulations["Basic mortality"].RunEcologicalProcess(gridCellCohorts, gridCellStocks, actingCohort,
//                cellEnvironment, deltas, madingleyCohortDefinitions, madingleyStockDefinitions, currentTimestep, trackProcesses, ref partial,
//                specificLocations, outputDetail, currentMonth, initialisation);
//        }
//
/** \brief
Update the properties of the acting cohort and of the environmental biomass pools after running the ecological processes for a cohort
*/
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param actingCohort The acting cohort 
@param cellEnvironment The environment of the current grid cell 
@param deltas The sorted list of deltas for the current grid cell 
@param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
@param madingleyStockDefinitions The definitions for stock functional groups in the model 
@param currentTimestep The current model time step 
@param tracker A process tracker 
//         void UpdateEcology(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, int[] actingCohort, 
//            map<string, double[]> cellEnvironment, Dictionary<string, Dictionary<string, double>> deltas, FunctionalGroupDefinitions 
//            madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions, unsigned currentTimestep, ProcessTracker tracker)
//        {
//            // Apply the results of within-cell ecological processes
//            ApplyEcologicalProcessResults.UpdateAllEcology(gridCellCohorts, actingCohort, cellEnvironment, deltas, currentTimestep, tracker);
//
//        }
//    }
//}
#endif
