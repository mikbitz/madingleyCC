#ifndef METABOLISM_H
#define METABOLISM_H
#include <IEcologicalProcessWithinGridCells.h>
#include <IMetabolismImplementation.h>
#include <TMetabolismHeterotroph.h>
#include <TMetabolismEndotherm.h>
#include <TMetabolismEctotherm.h>
#include <ThreadLocked.h>
/** \file Metabolism.h
 * \brief the Metabolism header file
 */

//
//namespace Madingley
//{
/** \brief  Performs metabolism */
class Metabolism : public IEcologicalProcessWithinGridCell
    {
    public:
/** \brief The available implementations of the metabolism process*/
          map<string, IMetabolismImplementation*> Implementations;
//
/** \brief Constructor Metabolism: fills the list of available implementations of metabolism*/
        Metabolism(string globalModelTimeStepUnit)
       {
           // Add the basic endotherm metabolism implementation to the list of implementations
           MetabolismEndotherm* MetabolismEndothermImplementation = new MetabolismEndotherm(globalModelTimeStepUnit);
           Implementations["basic endotherm"]=MetabolismEndothermImplementation;

           // Add the basic ectotherm metabolism implementation to the list of implementations
           MetabolismEctotherm* MetabolismEctothermImplementation = new MetabolismEctotherm(globalModelTimeStepUnit);
           Implementations["basic ectotherm"]=MetabolismEctothermImplementation;
       }
    ~Metabolism() {
        delete Implementations["basic endotherm"];
        delete Implementations["basic ectotherm"];
    }
/** \brief Initializes an implementation of metabolism
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
@param madingleyStockDefinitions The definitions for stock  functional groups in the model 
@param implementationKey The name of the implementation of metabolism to initialize  */
        void InitializeEcologicalProcess(GridCellCohortHandler& gridCellCohorts, GridCellStockHandler& gridCellStocks,
           FunctionalGroupDefinitions& madingleyCohortDefinitions, FunctionalGroupDefinitions& madingleyStockDefinitions,
           string implementationKey)
       {

       }

/** \brief Run metabolism
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
@param cellEnvironment The environment in the current grid cell 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
@param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
@param madingleyStockDefinitions The definitions for stock functional groups in the model 
@param currentTimestep The current model time step 
@param trackProcesses An instance of ProcessTracker to hold diagnostics for metabolism 
@param partial Thread-locked variables 
@param specificLocations Whether the model is being run for specific locations 
@param outputDetail The level of output detail being used for the current model run 
@param currentMonth The current model month  */
        void RunEcologicalProcess(GridCellCohortHandler& gridCellCohorts, GridCellStockHandler& gridCellStocks, 
           vector<int>& actingCohort, map<string, vector<double> >& cellEnvironment, map<string, map<string, double>>& deltas, 
           FunctionalGroupDefinitions& madingleyCohortDefinitions, FunctionalGroupDefinitions& madingleyStockDefinitions, 
           unsigned currentTimestep, ProcessTracker& trackProcesses, ThreadLockedParallelVariables& partial,
           bool specificLocations, string outputDetail, unsigned currentMonth, MadingleyModelInitialisation& initialisation)
       {
           double Realm = cellEnvironment["Realm"][0];
           if (madingleyCohortDefinitions.GetTraitNames("Heterotroph/Autotroph", gridCellCohorts[actingCohort].FunctionalGroupIndex) == "heterotroph")
           {
               if (madingleyCohortDefinitions.GetTraitNames("Endo/Ectotherm", gridCellCohorts[actingCohort].FunctionalGroupIndex) == "endotherm")
               {

                       Implementations["basic endotherm"]->RunMetabolism(gridCellCohorts, gridCellStocks, actingCohort, cellEnvironment, deltas, madingleyCohortDefinitions, madingleyStockDefinitions, currentTimestep, currentMonth);
               }
               else
               {
                       Implementations["basic ectotherm"]->RunMetabolism(gridCellCohorts, gridCellStocks, actingCohort, cellEnvironment, deltas, madingleyCohortDefinitions, madingleyStockDefinitions, currentTimestep, currentMonth);

               }

           }

           // If the process tracker is on and output detail is set to high and this cohort has not been merged yet, then record
           // the number of individuals that have died
 /*          if (trackProcesses.TrackProcesses && (outputDetail == "high"))
           {

               trackProcesses.TrackTimestepMetabolism((unsigned)cellEnvironment["LatIndex"][0],
                                               (unsigned)cellEnvironment["LonIndex"][0],
                                               currentTimestep,
                                               gridCellCohorts[actingCohort].IndividualBodyMass,
                                               actingCohort[0],
                                               cellEnvironment["Temperature"][currentMonth],
                                               deltas["biomass"]["metabolism"]);

           }
*/       }


    };
//}
#endif
