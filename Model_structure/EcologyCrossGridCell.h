#ifndef ECOLOGYCROSSGRID_H
#define ECOLOGYCROSSGRID_H
#include <IEcologicalProcessAcrossGridCells.h>
#include <IDispersalImplementation.h>
#include <ApplyCrossGridCellEcology.h>
#include <Dispersal.h>
//
//namespace Madingley
//{
/** \brief A class to specify, initalise and run ecological processes across grid cells */
class EcologyCrossGridCell
    {
    public:

/** \brief A vector of stopwatch objects for timing the ecological processes */
       vector<StopWatch> s2;              
/** \brief A sorted list of formulations of dispersal */
       map<string, IEcologicalProcessAcrossGridCells*>  DispersalFormulations;
/** \brief Get the sorted list of dispersal formulations */

/** \brief An instance of apply cross grid cell ecology */
        ApplyCrossGridCellEcology ApplyCrossGridCellEcologicalProcessResults;

/** \brief
Initalise the ecological processes
*/
void InitializeCrossGridCellEcology(string globalModelTimeStepUnit, bool drawRandomly, MadingleyModelInitialisation& modelInitialisation)
        {
           // Initialise dispersal formulations
           //_DispersalFormulations = new SortedList<string, IEcologicalProcessAcrossGridCells>();

           // Declare and attach dispersal formulations
           Dispersal *DispersalFormulation = new Dispersal(drawRandomly, globalModelTimeStepUnit, modelInitialisation);
           DispersalFormulations["Basic dispersal"]=DispersalFormulation;

           // Initialise apply ecology
           //ApplyCrossGridCellEcologicalProcessResults = new ApplyCrossGridCellEcology();
    }

    ~EcologyCrossGridCell() {
        delete DispersalFormulations["Basic dispersal"];
    }
/** \brief
Run ecological processes that operate across grid cells, for a particular grid cell. These should always occur after the within grid cell processes
*/

void RunCrossGridCellEcology(vector<unsigned>& cellIndex, ModelGrid& gridForDispersal, bool dispersalOnly, FunctionalGroupDefinitions& madingleyCohortDefinitions, FunctionalGroupDefinitions& madingleyStockDefinitions, uint currentMonth)
       {
           // RUN DISPERSAL
           DispersalFormulations["Basic dispersal"]->RunCrossGridCellEcologicalProcess(cellIndex, gridForDispersal, dispersalOnly, madingleyCohortDefinitions, madingleyStockDefinitions, currentMonth);       
       }

/** \brief
Update the properties of all cohorts across all grid cells
*/

void UpdateCrossGridCellEcology(ModelGrid& madingleyModelGrid, unsigned& dispersalCounter)
       {
           // Apply the results of cross-cell ecological processes
           ApplyCrossGridCellEcologicalProcessResults.UpdateAllCrossGridCellEcology(madingleyModelGrid, dispersalCounter);

       }
//   }
};
#endif