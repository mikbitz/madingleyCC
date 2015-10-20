#ifndef DISPERSAL_H
#define DISPERSAL_H
#include <string>
#include <IDispersalImplementation.h>
#include <GridCellCohortHandler.h>
#include <limits>
#include <TAdvectiveDispersal.h>
#include <TResponsiveDispersal.h>
#include <TDiffusiveDispersal.h>
using namespace std;
/** \file Dispersal.h
 * \brief the Dispersal header file
 */

/** \brief Performs dispersal */
class Dispersal  {
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    vector<std::reference_wrapper<Cohort>>disperseMonkeys;
    /** \brief The available implementations of the dispersal process */
    AdvectiveDispersal* AdvectiveDispersalImplementation;
    DiffusiveDispersal* DiffusiveDispersalImplementation;
    ResponsiveDispersal* ResponsiveDispersalImplementation;
    /** \brief Threshold (g) below which a marine individual is considered to be planktonic (i.e. cannot swim against the currents). Currently set to 10mg. */
    double PlanktonThreshold;
    public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
        Dispersal(){;}
    //----------------------------------------------------------------------------------------------
    /** \brief Setup for dispersal assigns pointers to current possible dispersal methods*/
    void setup(bool DrawRandomly, string globalModelTimeStepUnit, double PlanktonDispersalThreshold) {

        // Assign advective dispersal implementation
        AdvectiveDispersalImplementation = new AdvectiveDispersal(globalModelTimeStepUnit, DrawRandomly);

        // Assign diffusive dispersal implementation 
        DiffusiveDispersalImplementation = new DiffusiveDispersal(globalModelTimeStepUnit, DrawRandomly);

        // Assign responsive dispersal implementation 
        ResponsiveDispersalImplementation = new ResponsiveDispersal(globalModelTimeStepUnit, DrawRandomly);

        // Get the weight threshold below which organisms are dispersed planktonically
        PlanktonThreshold = PlanktonDispersalThreshold;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief tidy up pointers */
    ~Dispersal() {
    delete AdvectiveDispersalImplementation;
    delete DiffusiveDispersalImplementation;
    delete ResponsiveDispersalImplementation;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Run dispersal 
    @param cellIndex The cell index for the active cell in the model grid 
    @param gridForDispersal The model grid to run the process for 
    @param dispersalOnly Whether we are running dispersal only 
    @param madingleyCohortDefinitions The functional group definitions for cohorts in the model 
    @param madingleyStockDefinitions The functional group definitions for stocks in the model 
    @param currentMonth The current model month */
    void RunCrossGridCellEcologicalProcess(vector<unsigned>& cellIndex, ModelGrid& gridForDispersal, bool dispersalOnly, FunctionalGroupDefinitions& madingleyCohortDefinitions, FunctionalGroupDefinitions& madingleyStockDefinitions, unsigned currentMonth) {

        // Get the lat and lon indices
        unsigned ii = cellIndex[0];
        unsigned jj = cellIndex[1];

        // A boolean to check that the environmental layer exists
        bool varExists;
           
        // Check to see if the cell is marine
        double CellRealm = gridForDispersal.GetEnviroLayer("Realm", 0, ii, jj, varExists);

        // Go through all of the cohorts in turn and see if they disperse
        GridCellCohortHandler& WorkingGridCellCohorts = gridForDispersal.GetGridCellCohorts(ii, jj);
        // Loop through functional groups, and perform dispersal according to cohort type and status
        for (int kk = 0; kk < WorkingGridCellCohorts.size(); kk++) {
            // Work through the list of cohorts - be sure to work backward to get order of deletion right later
            for (int ll = WorkingGridCellCohorts[kk].size()-1;ll>0; ll--) {
                // Check to see if the cell is marine and the cohort type is planktonic
                bool dispersed=false;
                if (CellRealm == 2.0 &&
                        ((madingleyCohortDefinitions.GetTraitNames("Mobility", WorkingGridCellCohorts[kk][ll].FunctionalGroupIndex) == "planktonic") || (WorkingGridCellCohorts[kk][ll].IndividualBodyMass <= PlanktonThreshold))) {
                    // Run advective dispersal
                    AdvectiveDispersalImplementation->RunDispersal(disperseMonkeys,cellIndex, gridForDispersal, WorkingGridCellCohorts[kk][ll], kk, ll, currentMonth);
                }   // Otherwise, if mature do responsive dispersal
                else if (WorkingGridCellCohorts[kk][ll].MaturityTimeStep < std::numeric_limits<unsigned>::max()) {
                     //Run responsive dispersal
                    ResponsiveDispersalImplementation->RunDispersal(disperseMonkeys, cellIndex, gridForDispersal, WorkingGridCellCohorts[kk][ll], kk, ll, currentMonth);
                }    // If the cohort is immature, run diffusive dispersal
                else {
                    DiffusiveDispersalImplementation->RunDispersal(disperseMonkeys, cellIndex, gridForDispersal, WorkingGridCellCohorts[kk][ll], kk, ll, currentMonth);
                }
            }
        }
        // IF THE CELL IS MARINE, RUN ADVECTIVE DISPERSAL FOR THE PHYTOPLANKTON STOCK AS WELL IN v1
    }
    //----------------------------------------------------------------------------------------------
    void UpdateCrossGridCellEcology(ModelGrid& gridForDispersal, unsigned& dispersalCounter){
        dispersalCounter = disperseMonkeys.size();
        //do removals first as this currently depends closely on the order of cohorts in the grid
        for (auto& c: disperseMonkeys){
         gridForDispersal.DeleteGridCellIndividualCohort(c);
        }
        for (auto& c: disperseMonkeys){
         gridForDispersal.AddNewCohortToGridCell(c);
        }
        disperseMonkeys.clear();
    }
};

#endif
