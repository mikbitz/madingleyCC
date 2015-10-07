#ifndef STOCK_H
#define STOCK_H
/** \file Stock.h
 * \brief the Stock header file
 */

/** \brief
//    /// Hold individual stocks */
class Stock {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------

    /** \brief The index of the functional group that the stock belongs to */
    unsigned char FunctionalGroupIndex;

    /** \brief The mean body mass of an individual in this stock */

    double IndividualBodyMass;

    /** \brief The total biomass of the stock */
    double TotalBiomass;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    /** \brief default constructor - used only at definition in continain classes*/
    Stock() {
        ;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Constructor for stock class. Assigns stock starting properties
    @param functionalGroupIndex The functional group index of the stock being generated 
    @param individualMass The individual mass of the stock 
    @param initialTotalBiomass The initial total biomass of the stock */
    Stock(unsigned char functionalGroupIndex, double individualMass, double initialTotalBiomass) {
        FunctionalGroupIndex = functionalGroupIndex;
        IndividualBodyMass = individualMass;
        TotalBiomass = initialTotalBiomass;
    }
    //----------------------------------------------------------------------------------------------

};
#endif
