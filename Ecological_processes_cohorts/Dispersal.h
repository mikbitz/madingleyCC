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
    vector<Cohort>disperseMonkeys;
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
    void setup(MadingleyModelInitialisation& params) {

        // Assign advective dispersal implementation
        AdvectiveDispersalImplementation = new AdvectiveDispersal(params.GlobalModelTimeStepUnit, params.DrawRandomly);

        // Assign diffusive dispersal implementation 
        DiffusiveDispersalImplementation = new DiffusiveDispersal(params.GlobalModelTimeStepUnit, params.DrawRandomly);

        // Assign responsive dispersal implementation 
        ResponsiveDispersalImplementation = new ResponsiveDispersal(params.GlobalModelTimeStepUnit, params.DrawRandomly);

        // Get the weight threshold below which organisms are dispersed planktonically
        PlanktonThreshold = params.PlanktonDispersalThreshold;
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
    void RunCrossGridCellEcologicalProcess(GridCell& gcl, ModelGrid& gridForDispersal,  MadingleyModelInitialisation& params,  unsigned currentMonth) {

        gcl.ask([&](Cohort& c){
                // Check to see if the cell is marine and the cohort type is planktonic
                bool dispersed=false;
                
                if (gcl.isMarine() &&
                        ((params.CohortFunctionalGroupDefinitions.GetTraitNames("Mobility", c.FunctionalGroupIndex) == "planktonic") || (c.IndividualBodyMass <= PlanktonThreshold))) {
                    // Run advective dispersal
                    AdvectiveDispersalImplementation->RunDispersal(disperseMonkeys, gridForDispersal, c, currentMonth);
                }   // Otherwise, if mature do responsive dispersal
                else if (c.isMature()) {
                     //Run responsive dispersal
                    ResponsiveDispersalImplementation->RunDispersal(disperseMonkeys, gridForDispersal, c, currentMonth);
                }    // If the cohort is immature, run diffusive dispersal
                else {
                    DiffusiveDispersalImplementation->RunDispersal(disperseMonkeys, gridForDispersal, c,  currentMonth);
                }
            
        });
        // IF THE CELL IS MARINE, RUN ADVECTIVE DISPERSAL FOR THE PHYTOPLANKTON STOCK AS WELL IN v1
    }
    //----------------------------------------------------------------------------------------------
    void UpdateCrossGridCellEcology(ModelGrid& gridForDispersal, unsigned& dispersalCounter){
        dispersalCounter = disperseMonkeys.size();
        //do removals first as this currently depends closely on the order of cohorts in the grid
        //dispersals have been created in increasing order of index, so first reverse
        reverse(disperseMonkeys.begin(),disperseMonkeys.end());
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
