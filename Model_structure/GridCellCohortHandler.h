#ifndef GRIDCELLCOHORTHANDLER_H
#define GRIDCELLCOHORTHANDLER_H
#include <Cohort.h>
/** \file GridCellCohortHandler.h
 * \brief the GridCellCohortHandler header file
 */

/** \brief Handles the cohorts in a grid cell */

class GridCellCohortHandler
{
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    /** \brief A list of cohorts in the grid cell */
    vector< vector<Cohort> > GridCellCohorts;
public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    vector<Cohort>& operator[](unsigned i) {
        return GridCellCohorts[i];
    }
    //----------------------------------------------------------------------------------------------
    Cohort& operator[](vector<int> cll) {
        return GridCellCohorts[cll[0]][cll[1]];
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Create a new list of cohorts for the grid cell */
    GridCellCohortHandler() {
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Create a new list of cohorts of specified length corresponding to the number of functional groups 
    @param NumFunctionalGroups The number of functional groups for which there will be cohorts in this grid cell */
    GridCellCohortHandler(int NumFunctionalGroups){
        GridCellCohorts.resize(NumFunctionalGroups);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief If not setting the size in the constructor, this can be used instead
    @param NumFunctionalGroups The number of functional groups for which there will be cohorts in this grid cell */
    void setSize(int NumFunctionalGroups) {
        GridCellCohorts.resize(NumFunctionalGroups);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Gets the number of functional groups in this grid cell */
    unsigned size() {
        return GridCellCohorts.size();
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Gets the number of cohorts in this grid cell */
    int GetNumberOfCohorts() {
        int sum = 0;
        for (int ii = 0; ii < GridCellCohorts.size(); ii++) {
            sum += GridCellCohorts[ii].size();
        }
        return sum;
    }
    //----------------------------------------------------------------------------------------------
};
#endif
