#ifndef APPLYCROSSGRIDCELLECOLOGY_H
#define APPLYCROSSGRIDCELLECOLOGY_H
/** \file ApplyCrossGridCellEcology.h
 * \brief the ApplyCrossGridCellEcology header file
 */


//
//namespace Madingley
//{
/** \brief
Class for applying changes from the cross-grid cell ecological processes. These are held in matrices of lists in the modelgrid structure.
We simply loop through each cell, and check to see if there are any cohorts flagged as needing to be dispersed. If so, we point the cohort list
in the new grid cell to this cohort, we delete the pointer to it and to the new grid cell in the model grid delta structures, and we delete the pointer 
to it in the original cell

We can also output diagnostics here (temporariliy) as the whole grid needs to be completed before dispersal is enacted.*/
    class ApplyCrossGridCellEcology
   {
   public:
/** \brief
Apply all updates from the ecological processes to the properties of the acting cohort and to the environment
*/
        void UpdateAllCrossGridCellEcology(ModelGrid& madingleyModelGrid, unsigned& dispersalCounter)
       {
            dispersalCounter=0;
           // Loop through the delta array that holds the locations of the cohorts that are flagged as needing to be moved
           for (unsigned ii = 0; ii < madingleyModelGrid.DeltaFunctionalGroupDispersalArray.size() ; ii++)
           {
               for (unsigned jj = 0; jj < madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii].size(); jj++)
               {
                   // No cohorts to move if there are none in the delta dispersal array
                   if (madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii][ jj].size() == 0)
                   {
                   }
                   // Otherwise, loop through the cohorts and change the pointers/references to them one-by-one
                   else
                   {
                       for (int kk = 0; kk < madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii][ jj].size(); kk++)
                       {
                           // Find out which grid cell it is going to
                           vector<unsigned> CellToDisperseTo = madingleyModelGrid.DeltaCellToDisperseToArray[ii][ jj][kk];
                           
                           // Functional group is identified by the first array
                           unsigned CohortToDisperseFG = madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii][ jj][kk];

                           // Cohort number is identified by the second array
                           unsigned CohortToDisperseNum = madingleyModelGrid.DeltaCohortNumberDispersalArray[ii][ jj][kk];

                           // Simmply add it to the existing cohorts in that FG in the grid cell to disperse to
                           madingleyModelGrid.AddNewCohortToGridCell(CellToDisperseTo[0], CellToDisperseTo[1], (int)CohortToDisperseFG, madingleyModelGrid.GetGridCellIndividualCohort(ii, jj, (int)CohortToDisperseFG, (int)CohortToDisperseNum));

                           // Update the dispersal counter
                           dispersalCounter++;

                           // So now there is a pointer in the grid cell to which it is going. We have to delete the pointers in the original cell and in the
                           // delta array, but we need to do this without messing with the list structure; i.e. wait until all cohorts have been moved
                       }
                   }
               }
           }


           // Reset the delta array and remove the pointer to the cohort in the original list
           for (unsigned ii = 0; ii < madingleyModelGrid.DeltaFunctionalGroupDispersalArray.size(); ii++)
           {
               for (unsigned jj = 0; jj < madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii].size(); jj++)
               {
                   // No cohorts to move if there are none in the delta dispersal array
                   if (madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii][ jj].size() == 0)
                   {
                   }
                   // Otherwise, loop through the cohorts and change the pointers/references to them one-by-one
                   else
                   {
                       // Delete the cohorts from the original grid cell. Note that this needs to be done carefully to ensure that the correct ones 
                       // are deleted (lists shift about when an internal element is deleted.
                       madingleyModelGrid.DeleteGridCellIndividualCohorts(ii, jj, madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii][ jj], madingleyModelGrid.DeltaCohortNumberDispersalArray[ii][jj]);

                        // Reset the lists in the delta dispersal arrays
                       madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii][ jj].clear();
                       madingleyModelGrid.DeltaCohortNumberDispersalArray[ii][ jj].clear();

                       // Reset the list in the grid cells to disperse to array
                       madingleyModelGrid.DeltaCellToDisperseToArray[ii][ jj].clear();
                   }
               }
           }

       }
   };
//}
#endif
