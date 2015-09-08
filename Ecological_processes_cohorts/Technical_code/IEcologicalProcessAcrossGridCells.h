#ifndef IECOLOGICALPROCESSACROSSGRIDCELLS_H
#define IECOLOGICALPROCESSACROSSGRIDCELLS_H
#include <vector>
#include <ModelGrid.h>
#include <FunctionalGroupDefinitions.h>
/** \file IEcologicalProcessAcrossGridCells.h
 * \brief the IEcologicalProcessAcrossGridCells header file
 */
//C++ pure virtual class - acts as an interface
//namespace Madingley
//{
/** \brief Interface for cross grid-cell ecological process code */
class IEcologicalProcessAcrossGridCells
    {
    public:
/** \brief Run the cross-grid-cell ecological process 
@param cellIndex The cell index for the active cell in the model grid 
@param gridForDispersal The model grid to run the process for 
@param dispersalOnly Whether we are running dispersal only 
@param madingleyCohortDefinitions The functional group definitions for cohorts in the model 
@param madingleyStockDefinitions The functional group definitions for stocks in the model 
@param currentMonth The current model month */
virtual void RunCrossGridCellEcologicalProcess(vector<unsigned> cellIndex, ModelGrid gridForDispersal, bool dispersalOnly, 
           FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions, 
           unsigned currentMonth){;}
   };
//}
#endif
