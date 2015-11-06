#ifndef EATING_H
#define EATING_H
#include <IEatingImplementation.h>
#include <IEcologicalProcessWithinGridCells.h>
#include <random>
#include <TRevisedPredation.h>
#include <TRevisedHerbivory.h>
/** \file Eating.h
 * \brief the Eating header file
 */

/** \brief Performs eating */
class Eating : public IEcologicalProcessWithinGridCell {
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    /** \brief The available implementations of the eating process */
    map<string, IEatingImplementation*> Implementations;
    /** \brief Tracks the total time to handle all potential food for omnivores */
    double TotalTimeToEatForOmnivores;
    /** \brief An instance of the simple random number generator class */
    std::default_random_engine RandomNumberGenerator;
public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    /** \brief Constructor for Eating: fills the list of available implementations of eating */
    Eating(double cellArea, string globalModelTimeStepUnit) {

        // Add the revised herbivory implementation to the list of implementations
        RevisedHerbivory *RevisedHerbivoryImplementation = new RevisedHerbivory(cellArea, globalModelTimeStepUnit);
        Implementations["revised herbivory"] = RevisedHerbivoryImplementation;
        //Add the revised predation implementation to the list of implementations
        RevisedPredation *RevisedPredationImplementation = new RevisedPredation(cellArea, globalModelTimeStepUnit);
        Implementations["revised predation"] = RevisedPredationImplementation;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Destructor to remove pointer storage */
    ~Eating() {
        delete Implementations["revised herbivory"];
        delete Implementations["revised predation"];
    }
    //----------------------------------------------------------------------------------------------
    /** \briefInitializes an implementation of eating 
    @param gridCellCohorts The cohorts in the current grid cell 
    @param gridCellStocks The stocks in the current grid cell 
    @param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
    @param madingleyStockDefinitions The definitions for stock functional groups in the model 
    @param implementationKey The name of the implementation of eating to initialize 
    \remarks Eating needs to be initialized every time step */
    void InitializeEcologicalProcess(GridCellCohortHandler& gridCellCohorts, GridCellStockHandler& gridCellStocks,
            FunctionalGroupDefinitions& madingleyCohortDefinitions, FunctionalGroupDefinitions& madingleyStockDefinitions,
            string implementationKey) {

        // Initialize the implementation of the eating process
        Implementations[implementationKey]->InitializeEatingPerTimeStep(gridCellCohorts, gridCellStocks,
                madingleyCohortDefinitions, madingleyStockDefinitions);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Run eating 
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
    @param outputDetail The level of output detail being used for the current model run 
    @param currentMonth The current model month  */
    void RunEcologicalProcess(GridCellCohortHandler& gridCellCohorts,
            GridCellStockHandler& gridCellStocks, Cohort& actingCohort,
            map<string, vector<double> >& cellEnvironment,
            map<string, map<string, double> >& deltas,
            FunctionalGroupDefinitions& madingleyCohortDefinitions,
            FunctionalGroupDefinitions& madingleyStockDefinitions,
            unsigned currentTimestep, ProcessTracker& trackProcesses,
            ThreadLockedParallelVariables& partial,
            unsigned currentMonth, MadingleyModelInitialisation& initialisation) {
        // Get the nutrition source (herbivory, carnivory or omnivory) of the acting cohort
        string NutritionSource = madingleyCohortDefinitions.GetTraitNames("Nutrition source", actingCohort.FunctionalGroupIndex);
        map<string, int> vores;
        vores["herbivore"] = 0;
        vores["carnivore"] = 1;
        vores["omnivore" ] = 2;

        //            // Switch to the appropriate eating process(es) given the cohort's nutrition source
        switch (vores[NutritionSource]) {
            case 0://"herbivore":

                // Get the assimilation efficiency for herbivory for this cohort from the functional group definitions
                Implementations["revised herbivory"]->AssimilationEfficiency =
                        madingleyCohortDefinitions.GetBiologicalPropertyOneFunctionalGroup
                        ("herbivory assimilation", actingCohort.FunctionalGroupIndex);

                // Get the proportion of time spent eating for this cohort from the functional group definitions
                Implementations["revised herbivory"]->ProportionTimeEating = actingCohort.ProportionTimeActive;

                // Calculate the potential biomass available from herbivory
                if (cellEnvironment["Realm"][0] == 2.0)
                    Implementations["revised herbivory"]->GetEatingPotentialMarine
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, madingleyCohortDefinitions, madingleyStockDefinitions);
                else

                    Implementations["revised herbivory"]->GetEatingPotentialTerrestrial
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, madingleyCohortDefinitions, madingleyStockDefinitions);

                // Run herbivory to apply changes in autotroph biomass from herbivory and add biomass eaten to the delta arrays
                Implementations["revised herbivory"]->RunEating
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, deltas, madingleyCohortDefinitions,
                        madingleyStockDefinitions, trackProcesses,
                        currentTimestep,  initialisation);

                break;

            case 1://"carnivore":

                // Get the assimilation efficiency for predation for this cohort from the functional group definitions
                Implementations["revised predation"]->AssimilationEfficiency =
                        madingleyCohortDefinitions.GetBiologicalPropertyOneFunctionalGroup
                        ("carnivory assimilation", actingCohort.FunctionalGroupIndex);

                Implementations["revised predation"]->ProportionTimeEating = actingCohort.ProportionTimeActive;

                // Calculate the potential biomass available from predation
                if (cellEnvironment["Realm"][0] == 2.0)
                    Implementations["revised predation"]->GetEatingPotentialMarine
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, madingleyCohortDefinitions, madingleyStockDefinitions);
                else
                    Implementations["revised predation"]->GetEatingPotentialTerrestrial
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, madingleyCohortDefinitions, madingleyStockDefinitions);
                // Run predation to apply changes in prey biomass from predation and add biomass eaten to the delta arrays
                Implementations["revised predation"]->RunEating
                        (gridCellCohorts, gridCellStocks, actingCohort, cellEnvironment, deltas,
                        madingleyCohortDefinitions, madingleyStockDefinitions, trackProcesses,
                        currentTimestep,  initialisation);


                break;

            case 2://"omnivore":

                // Get the assimilation efficiency for predation for this cohort from the functional group definitions 
                Implementations["revised predation"]->AssimilationEfficiency =
                        madingleyCohortDefinitions.GetBiologicalPropertyOneFunctionalGroup
                        ("carnivory assimilation", actingCohort.FunctionalGroupIndex);

                // Get the assimilation efficiency for herbivory for this cohort from the functional group definitions
                Implementations["revised herbivory"]->AssimilationEfficiency =
                        madingleyCohortDefinitions.GetBiologicalPropertyOneFunctionalGroup
                        ("herbivory assimilation", actingCohort.FunctionalGroupIndex);

                // Get the proportion of time spent eating and assign to both the herbivory and predation implementations
                //                   double ProportionTimeEating = actingCohort.ProportionTimeActive;
                Implementations["revised predation"]->ProportionTimeEating = actingCohort.ProportionTimeActive;
                Implementations["revised herbivory"]->ProportionTimeEating = actingCohort.ProportionTimeActive;

                // Calculate the potential biomass available from herbivory
                if (cellEnvironment["Realm"][0] == 2.0)
                    Implementations["revised herbivory"]->GetEatingPotentialMarine
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, madingleyCohortDefinitions,
                        madingleyStockDefinitions);
                else
                    Implementations["revised herbivory"]->GetEatingPotentialTerrestrial
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, madingleyCohortDefinitions,
                        madingleyStockDefinitions);

                // Calculate the potential biomass available from predation
                if (cellEnvironment["Realm"][0] == 2.0)
                    Implementations["revised predation"]->GetEatingPotentialMarine
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, madingleyCohortDefinitions,
                        madingleyStockDefinitions);
                else
                    Implementations["revised predation"]->GetEatingPotentialTerrestrial
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, madingleyCohortDefinitions,
                        madingleyStockDefinitions);

                // Calculate the total handling time for all expected kills from predation and expected plant matter eaten in herbivory
                TotalTimeToEatForOmnivores =
                        Implementations["revised herbivory"]->TimeUnitsToHandlePotentialFoodItems +
                        Implementations["revised predation"]->TimeUnitsToHandlePotentialFoodItems;

                // Assign this total time to the relevant variables in both herbviory and predation, so that actual amounts eaten are calculated correctly
                Implementations["revised herbivory"]->TimeUnitsToHandlePotentialFoodItems = TotalTimeToEatForOmnivores;
                Implementations["revised predation"]->TimeUnitsToHandlePotentialFoodItems = TotalTimeToEatForOmnivores;

                // Run predation to update prey cohorts and delta biomasses for the acting cohort
                Implementations["revised predation"]->RunEating
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, deltas, madingleyCohortDefinitions,
                        madingleyStockDefinitions, trackProcesses,
                        currentTimestep, initialisation);

                // Run herbivory to update autotroph biomass and delta biomasses for the acting cohort
                Implementations["revised herbivory"]->RunEating
                        (gridCellCohorts, gridCellStocks, actingCohort,
                        cellEnvironment, deltas, madingleyCohortDefinitions,
                        madingleyStockDefinitions, trackProcesses,
                        currentTimestep, initialisation);

                break;

            default:

                // For nutrition source that are not supported, throw an error
                cout << "The model currently does not contain an eating model for nutrition source:" << NutritionSource << endl;
                exit(1);
                break;

        }
        // Check that the biomasses from predation and herbivory in the deltas is a number
        assert(!std::isnan(deltas["biomass"]["predation"]) && "BiomassFromEating is NaN");
        assert(!std::isnan(deltas["biomass"]["herbivory"]) && "BiomassFromEating is NaN");

    }
    //----------------------------------------------------------------------------------------------
};
#endif
