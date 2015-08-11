#ifndef ECOLOGYCROSSGRID_H
#define ECOLOGYCROSSGRID_H
//
//namespace Madingley
//{
/** \brief A class to specify, initalise and run ecological processes across grid cells */
class EcologyCrossGridCell
    {
    public:
//        # region Properties and Fields
//
/** \brief
A vector of stopwatch objects for timing the ecological processes
*/
//        public StopWatch[] s2;
//               
/** \brief
A sorted list of formulations of dispersal
*/
//        private SortedList<string, IEcologicalProcessAcrossGridCells> _DispersalFormulations;
/** \brief
Get the sorted list of dispersal formulations
*/
//        public SortedList<string, IEcologicalProcessAcrossGridCells> DispersalFormulations
//	    {
//            get { return _DispersalFormulations; }
//        }
//	
/** \brief
An instance of apply cross grid cell ecology
*/
//        ApplyCrossGridCellEcology ApplyCrossGridCellEcologicalProcessResults;
//
//
//        # endregion
//
/** \brief
Initalise the ecological processes
*/
void InitializeCrossGridCellEcology(string globalModelTimeStepUnit, bool drawRandomly, MadingleyModelInitialisation modelInitialisation)
        {
//            // Initialise dispersal formulations
//            _DispersalFormulations = new SortedList<string, IEcologicalProcessAcrossGridCells>();
//
//            // Declare and attach dispersal formulations
//            Dispersal DispersalFormulation = new Dispersal(drawRandomly, globalModelTimeStepUnit, modelInitialisation);
//            _DispersalFormulations.Add("Basic dispersal", DispersalFormulation);
//
//            // Initialise apply ecology
//            ApplyCrossGridCellEcologicalProcessResults = new ApplyCrossGridCellEcology();
        }

/** \brief
Run ecological processes that operate across grid cells, for a particular grid cell. These should always occur after the within grid cell processes
*/
//
//        public void RunCrossGridCellEcology(uint[] cellIndex, ModelGrid gridForDispersal, bool dispersalOnly, FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions, uint currentMonth)
//        {
//            // RUN DISPERSAL
//            _DispersalFormulations["Basic dispersal"].RunCrossGridCellEcologicalProcess(cellIndex, gridForDispersal, dispersalOnly, madingleyCohortDefinitions, madingleyStockDefinitions, currentMonth);       
//        }
//
/** \brief
Update the properties of all cohorts across all grid cells
*/
//
//        public void UpdateCrossGridCellEcology(ModelGrid madingleyModelGrid, ref uint dispersalCounter)
//        {
//            // Apply the results of cross-cell ecological processes
//            ApplyCrossGridCellEcologicalProcessResults.UpdateAllCrossGridCellEcology(madingleyModelGrid, ref dispersalCounter);
//
//        }
//    }
};
#endif