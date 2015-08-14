#ifndef STOCK_H
#define STOCK_H
/** \file Stock.h
 * \brief the Stock header file
 */


//
//namespace Madingley
//{
/** \brief
//    /// Hold individual stocks */
class Stock
    {
    public:
/** \brief
The index of the functional group that the stock belongs to
*/
     unsigned char FunctionalGroupIndex;
/** \brief
Get and set the functional group that the stock belongs to
*/
//         byte FunctionalGroupIndex { get { return _FunctionalGroupIndex; } }
//
/** \brief
The mean body mass of an individual in this stock
*/
//        private double _IndividualBodyMass;
/** \brief
Get and set the mean body mass of an individual in this stock
*/
         double IndividualBodyMass;
//        {
//            get { return _IndividualBodyMass; }
//            set { _IndividualBodyMass = value; }
//        }
//
/** \brief
The total biomass of the stock
*/
        double TotalBiomass;
/** \brief
Get and set the total biomass of this stock
*/
//         double TotalBiomass
//        {
//            get { return _TotalBiomass; }
//            set { _TotalBiomass = value; }
//        }
//
Stock(){;}
/** \brief
Constructor for stock class. Assigns stock starting properties
@param functionalGroupIndex The functional group index of the stock being generated 
@param individualMass The individual mass of the stock 
@param initialTotalBiomass The initial total biomass of the stock */
        Stock(unsigned char functionalGroupIndex, double individualMass, double initialTotalBiomass)
       {
           FunctionalGroupIndex = functionalGroupIndex;
           IndividualBodyMass = individualMass;
           TotalBiomass = initialTotalBiomass;
       }
       

    };
//}
#endif
