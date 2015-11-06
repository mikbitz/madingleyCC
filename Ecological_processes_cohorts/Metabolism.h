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
class Metabolism : public IEcologicalProcessWithinGridCell {
public:
    /** \brief The available implementations of the metabolism process*/
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    map<string, IMetabolismImplementation*> Implementations;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    /** \brief Constructor Metabolism: fills the list of available implementations of metabolism*/
    Metabolism(string globalModelTimeStepUnit) {
        // Add the basic endotherm metabolism implementation to the list of implementations
        MetabolismEndotherm* MetabolismEndothermImplementation = new MetabolismEndotherm(globalModelTimeStepUnit);
        Implementations["basic endotherm"] = MetabolismEndothermImplementation;

        // Add the basic ectotherm metabolism implementation to the list of implementations
        MetabolismEctotherm* MetabolismEctothermImplementation = new MetabolismEctotherm(globalModelTimeStepUnit);
        Implementations["basic ectotherm"] = MetabolismEctothermImplementation;
    }
    //----------------------------------------------------------------------------------------------
    /** Destrcutor to tidy up pointers*/
    ~Metabolism() {
        delete Implementations["basic endotherm"];
        delete Implementations["basic ectotherm"];
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Initializes an implementation of metabolism
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock  functional groups in the model 
    @param implementationKey The name of the implementation of metabolism to initialize  */
    void InitializeEcologicalProcess(GridCellCohortHandler& gridCellCohorts, GridCellStockHandler& gridCellStocks,
            FunctionalGroupDefinitions& madingleyCohortDefinitions, FunctionalGroupDefinitions& madingleyStockDefinitions,
            string implementationKey) {

    }
    //----------------------------------------------------------------------------------------------
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
            Cohort& actingCohort, map<string, vector<double> >& cellEnvironment, map<string, map<string, double>>&deltas,
            FunctionalGroupDefinitions& madingleyCohortDefinitions, FunctionalGroupDefinitions& madingleyStockDefinitions,
            unsigned currentTimestep, ThreadLockedParallelVariables& partial,
            unsigned currentMonth, MadingleyModelInitialisation& initialisation) {
        double Realm = cellEnvironment["Realm"][0];
        if (madingleyCohortDefinitions.GetTraitNames("Heterotroph/Autotroph", actingCohort.FunctionalGroupIndex) == "heterotroph") {
            if (madingleyCohortDefinitions.GetTraitNames("Endo/Ectotherm", actingCohort.FunctionalGroupIndex) == "endotherm") {

                Implementations["basic endotherm"]->RunMetabolism(  actingCohort, cellEnvironment, deltas, currentTimestep, currentMonth);
            } else {
                Implementations["basic ectotherm"]->RunMetabolism(  actingCohort, cellEnvironment, deltas, currentTimestep, currentMonth);

            }

        }

    }
    //----------------------------------------------------------------------------------------------
};
#endif
