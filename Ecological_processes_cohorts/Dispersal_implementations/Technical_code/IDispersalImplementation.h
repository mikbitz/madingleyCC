#ifndef IDISPERSALIMPLEMENTATION_H
#define IDISPERSALIMPLEMENTATION_H
#include <ModelGrid.h>
#include <Cohort.h>

/** \file IDispersalImplementation.h
 * \brief the IDispersalImplementation header file
 */

//
//namespace Madingley
//{
/** \brief Interface for implementations of the ecological process of dispersal */
//C++ virtual class acts as an interface - however this may mean we need pointers...
class IDispersalImplementation
    {
    public:
/** \brief Time units associated with the formulation of dispersal */
     //string TimeUnitImplementation;
//        
/** \brief Scalar to convert from time units associated with dispersal to the global model time step unit */
    //DoubleProperty DeltaT;
//
/** \brief Run the dispersal implementation */
     virtual void RunDispersal(vector<unsigned>& cellIndex, ModelGrid& gridForDispersal, Cohort& cohortToDisperse, int actingCohortFunctionalGroup, int actingCohortNumber, unsigned currentMonth){;}
    };

//}
#endif
