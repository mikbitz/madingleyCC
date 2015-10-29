#ifndef GRIDCELLSTOCKHANDLER_H
#define GRIDCELLSTOCKHANDLER_H
#include <Stock.h>

/** \file GridCellStockHandler.h
 * \brief the GridCellStockHandler header file
 */

/** \brief Handles the stocks in a grid cell */

class GridCellStockHandler {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    /** \brief A vector (with elements correpsonding to functional groups) of lists of stocks in the current grid cell*/
    vector< vector<Stock> > GridCellStocks;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    vector<Stock>& operator[](unsigned i) {
        return GridCellStocks[i];
    }
    //----------------------------------------------------------------------------------------------
    Stock& operator[](vector<int> s) {
        return GridCellStocks[s[0]][s[1]];
    }
    //----------------------------------------------------------------------------------------------
    unsigned size() {
        return GridCellStocks.size();
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Overloaded constructor for the grid cell stock handler: initialises a new vector of lists of stocks*/
    GridCellStockHandler() {}
    //----------------------------------------------------------------------------------------------
    /** \brief Overloaded constructor for the grid cell stock handler: initialises a new vector of lists of stocks with number of elements equal to the number of functional groups 
    @param NumFunctionalGroups The number of stock functional groups in the model */
    GridCellStockHandler(int NumFunctionalGroups){
        GridCellStocks.resize(NumFunctionalGroups);
    }
    //----------------------------------------------------------------------------------------------
    void setSize(int NumFunctionalGroups) {
        GridCellStocks.resize(NumFunctionalGroups);
    }
    //----------------------------------------------------------------------------------------------
 
};

#endif
