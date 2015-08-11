#ifndef IREPRODUCTIONIMPLEMENTATION_H
#define IREPRODUCTIONIMPLEMENTATION_H
#include <Properties.h>
#include <GridCellCohortHandler.h>
#include <GridCellStockHandler.h>
#include <FunctionalGroupDefinitions.h>
#include <ProcessTracker.h>
#include <MadingleyModelInitialisation.h>
#include <map>
/** \file IReproductionImplementation.h
 * \brief the IReproductionImplementation header file
 */

//
//namespace Madingley
//{
/** \brief Interface for implementations of the ecological process of reproduction */
class IReproductionImplementation
    {
    public:
//        #region Property and field definitions
//        
/** \brief
Time units associated with the formulation of reproduction
*/
      //StringProperty TimeUnitImplementation;

/** \brief
Scalar to convert from the time units associated with reproduction to the global model time step unit
*/
      //DoubleProperty DeltaT;
       
//        #endregion
//
/** \brief Generate new cohorts from reproductive potential mass
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
@param cellEnvironment The environment in the current grid cell 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
@param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
@param madingleyStockDefinitions The definitions for stock functional groups in the model 
@param currentTimestep The current model time step 
@param trackProcesses An instance of ProcessTracker to hold diagnostics for eating 
@param partial Thread-locked variables 
@param iteroparous Whether the acting cohort is iteroparous, as opposed to semelparous 
@param currentMonth The current model month */
       virtual void RunReproductionEvents(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, vector<int> actingCohort, 
           map<string, vector<double> > cellEnvironment, map<string, map<string, double> > deltas, FunctionalGroupDefinitions 
           madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions, unsigned currentTimestep, ProcessTracker trackProcesses, 
           ThreadLockedParallelVariables& partial, bool iteroparous, unsigned currentMonth);
       
/** \brief Assigns surplus body mass to reproductive potential mass
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
@param cellEnvironment The environment in the current grid cell 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
@param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
@param madingleyStockDefinitions The definitions for stock functional groups in the model 
@param currentTimestep The current model time step 
@param trackProcesses An instance of ProcessTracker to hold diagnostics for reproduction */
       virtual void RunReproductiveMassAssignment(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, vector<int> actingCohort, map<string, vector<double>> cellEnvironment, map<string, map<string, double> > deltas, FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions, unsigned currentTimestep, ProcessTracker trackProcesses);
    };
//}
#endif
