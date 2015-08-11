#ifndef IEATINGIMPLEMENTATION_H
#define IEATINGIMPLEMENTATION_H
#include <Properties.h>
#include <GridCellCohortHandler.h>
#include <GridCellStockHandler.h>
#include <FunctionalGroupDefinitions.h>
#include <ProcessTracker.h>
#include <MadingleyModelInitialisation.h>
#include <map>
/** \file IEatingImplementation.h
 * \brief the IEatingImplementation header file
 */


//namespace Madingley
//{
/** \brief //    /// Interface for implementations of the ecological process of eating */
class IEatingImplementation
    {
    public:
/** \brief Time units associated with the formulation of eating*/
        //StringProperty TimeUnitImplementation;       
/** \brief Scalar to convert from time units associated with eating to the global model time step unit*/
       // DoubleProperty DeltaT; 
/** \brief Assimilation efficiency of food mass into acting cohort mass*/
        double AssimilationEfficiency;
/** \brief Proportion of time spent eating*/
        double ProportionTimeEating;
/** \brief Time to handle all prey cohorts or plant mass encountered*/
        double TimeUnitsToHandlePotentialFoodItems;      
/** \brief List of functional group indices to act on*/
        vector<int> FunctionalGroupIndicesToEat;
/** \brief The total biomass eaten by the acting cohort */
        double TotalBiomassEatenByCohort;
/** \brief Initialises eating implementation each time step
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param madingleyCohortDefinitions The definitions for cohorts in the model 
@param madingleyStockDefinitions The definitions for stocks in the model  */ 
        void InitializeEatingPerTimeStep(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions);
//        
/** \brief Calculate the potential biomass that could be gained through eating for marine cells
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param actingCohort The position of the acting cohort in the jagged array of cohorts 
@param cellEnvironment The environment in the current grid cell 
@param madingleyCohortDefinitions The definitions for cohorts in the model 
@param madingleyStockDefinitions The definitions for stocks in the model  */
       virtual void GetEatingPotentialMarine(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, vector<int> actingCohort, 
           map<string, vector<double>> cellEnvironment, FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions 
           madingleyStockDefinitions);

/** \brief Calculate the potential biomass that could be gained through eating for terrestrial cells
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param actingCohort The position of the acting cohort in the jagged array of cohorts 
@param cellEnvironment The environment in the current grid cell 
@param madingleyCohortDefinitions The definitions for cohorts in the model 
@param madingleyStockDefinitions The definitions for stocks in the model  */
       virtual void GetEatingPotentialTerrestrial(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, vector<int> actingCohort,
           map<string, vector<double> > cellEnvironment, FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions
           madingleyStockDefinitions);


/** \brief Calculate the actual biomass eaten from each cohort or sotck, apply changes from eating to the cohorts or stocks eaten, and update deltas for the acting cohort
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
@param cellEnvironment The environment in the current grid cell 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
@param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
@param madingleyStockDefinitions The definitions for stock functional groups in the model 
@param trackProcesses An instance of ProcessTracker to hold diagnostics for eating 
@param currentTimestep The current model time step 
@param specificLocations Whether the model is being run for specific locations 
@param outputDetail The level of output detail being used in this model run  */ 
       virtual void RunEating(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, 
           vector<int> actingCohort, map<string, vector< double > > cellEnvironment, 
           map<string, map<string, double>> deltas, 
           FunctionalGroupDefinitions madingleyCohortDefinitions, 
           FunctionalGroupDefinitions madingleyStockDefinitions, 
           ProcessTracker trackProcesses, unsigned currentTimestep,
           bool specificLocations, string outputDetail, MadingleyModelInitialisation initialisation);

    };
//}
#endif