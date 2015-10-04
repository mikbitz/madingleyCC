#ifndef GRIDCELLCOHORTHANDLER_H
#define GRIDCELLCOHORTHANDLER_H
#include <Cohort.h>
/** \file GridCellCohortHandler.h
 * \brief the GridCellCohortHandler header file
 */

//
//namespace Madingley
//{
/** \brief Handles the cohorts in a grid cell

//    /// <todoD>NOTE TO DT: CAN USE COLLECTIONBASE SYNTAX TO ELIMINATE A LOT OF CODE HERE</todoD>
//    /// <todo>Create a wrapper class to handle our array of lists of gridCellCohorts within an individual grid cell</todo>
*/
//MB vector< vector<Cohort> > I think does the job?
class GridCellCohortHandler// : IList<List<Cohort>>, IEnumerable<List<Cohort>>
    {
/** \brief A list of cohorts in the grid cell */
vector< vector<Cohort> > GridCellCohorts;
    public:
vector<Cohort>& operator[](unsigned i){return GridCellCohorts[i];}
Cohort& operator[](vector<int> cll){return GridCellCohorts[cll[0]][cll[1]];}
unsigned size(){return GridCellCohorts.size();}
//        private List<Cohort> [] GridCellCohorts;
//
/** \brief Create a new list of cohorts for the grid cell */
        GridCellCohortHandler()
       {
//            GridCellCohorts = new List<Cohort>[0];
       }

/** \brief Create a new list of cohorts of specified length corresponding to the number of functional groups 
@param NumFunctionalGroups The number of functional groups for which there will be cohorts in this grid cell */
//        GridCellCohortHandler(int NumFunctionalGroups)
       void setSize(int NumFunctionalGroups)
       {
            GridCellCohorts.resize(NumFunctionalGroups);
       }

/** \brief Update grid cell cohorts with a specified list of cohorts 
@param ExistingCohorts A list of cohorts to update the grid cell cohorts with */
//         GridCellCohortHandler(List<Cohort>[] ExistingCohorts)
//        {
//            GridCellCohorts = ExistingCohorts; //MB !squeak! memory leak??
//        }
//
/** \brief Get or set the list of cohorts for a specified functional group index //MB well controlled?? leaky??
@param functionalGroupIndex The index of the functional group to get or set the list of cohorts for 
@return The list of cohorts in the specified functional group*/
//         List<Cohort> this[int functionalGroupIndex]
//        {
//            get { return GridCellCohorts[functionalGroupIndex]; }
//            set
//            {
//                GridCellCohorts[functionalGroupIndex] = value;
//            }
//        }
//
/** \brief Gets or sets a particular cohort within the grid cell cohorts
@param index A vector of two values corresponding to the functional group index and the index of the desired cohort within this functional group 
@return The specified cohort*/
//         Cohort this[int[] index]
//        {
//            get { return GridCellCohorts[index[0]][index[1]]; }
//            set { GridCellCohorts[index[0]][index[1]] = value; }
//        }
//
//         // Gets of sets a cohort within the array of lists of gridCellCohorts where the first element of the 2-element vector passed in is the array index and the second element is the list index
/** \brief Gets or sets a particular cohort within the grid cell cohorts
@param functionalGroupIndex The functional group index of the desired cohort 
@param cohortIndex The index of the cohort within the specified functional group 
@return The specified cohort*/
//         Cohort this[int functionalGroupIndex, int cohortIndex]
//        {
//            get { return GridCellCohorts[functionalGroupIndex][cohortIndex]; }
//            set
//            {
//                if (GridCellCohorts[functionalGroupIndex] == null) GridCellCohorts[functionalGroupIndex] = new List<Cohort>();
//                GridCellCohorts[functionalGroupIndex].Add(value);
//            }
//        }
//
//        
/** \brief Get the functional group index a specified cohort
@param cohort The cohort to return the functional group index for 
@return The functional group index of the specified cohort*/
//         int IndexOf(List<Cohort> cohort)
//        {
//            return ((IList<List<Cohort>>)GridCellCohorts).IndexOf(cohort);
//        }
//
/** \brief Inserts a new list of cohorts at a specified functional group index - CURRENTLY  NOT SUPPORTED

@param index The index in the list of functional groups to insert the list of cohorts in 
@param listOfCohorts The list of cohorts to insert */
//         void Insert(int index, List<Cohort> listOfCohorts)
//        {
//            Debug.Fail("The model does not currently support the addition of functional groups");
//            ((IList<List<Cohort>>)GridCellCohorts).Insert(index, listOfCohorts);
//        }
//
/** \brief Removes a list of cohorts in a specified functional group - CURRENTLY NOT SUPPORTED

@param functionalGroupIndex The index of the functional group to remove the list of cohorts for */
//         void RemoveAt(int functionalGroupIndex)
//        {
//            Debug.Fail("The model does not currently support the removal of functional groups");
//            ((IList<List<Cohort>>)GridCellCohorts).RemoveAt(functionalGroupIndex);
//        }
//
/** \brief Adds a list of cohorts at the end of the functional group indices - CURRENTLY NOT SUPPORTED
@param listOfCohorts The list of cohorts to add */
//         void Add(List<Cohort> listOfCohorts)
//        {
//            Debug.Fail("The model does not currently support the addition of functional groups");
//            ((IList<List<Cohort>>)GridCellCohorts).Add(listOfCohorts);
//        }
//
/** \brief Currently not implemented */
//         void Clear()
//        {
//            throw new NotImplementedException();
//        }
//
/** \brief Currently not implemented
@param item NA 
@return NA*/
//         bool Contains(List<Cohort> item)
//        {
//            throw new NotImplementedException();
//        }
//
/** \brief Currently not implemented
@param array NA 
@param arrayIndex NA 
*/
//         void CopyTo(List<Cohort>[] array, int arrayIndex)
//        {
//            throw new NotImplementedException();
//        }
//
/** \brief Gets the number of functional groups in the grid cell cohorts */
//         int Count
//        {
//            get { return GridCellCohorts.Count(); }
//        }
//
/** \brief Currently not implemented */
//         bool IsReadOnly
//        {
//            get { throw new NotImplementedException(); }
//        }
//
/** \brief Currently not implemented
@param item NA 
@return NA*/
//         bool Remove(List<Cohort> item)
//        {
//            throw new NotImplementedException();
//        }
//
/** \brief Returns an the grid cell cohorts as an IEnumerator 
@return The grid cell cohorts as an IEnumerator*/
//         IEnumerator<List<Cohort>> GetEnumerator()
//        {
//            return new GridCellCohortsEnum(GridCellCohorts);
//        }
//
/** \brief Return an IEnumerable as an IEnumerator
@return The IEnumerable as an IEnumerator*/
//        IEnumerator IEnumerable.GetEnumerator()
//        {
//            return (IEnumerator)GetEnumerator();
//        }
//
/** \brief Gets the number of cohorts in this grid cell */
         int GetNumberOfCohorts()
        {

            int sum = 0;
            for (int ii = 0; ii < GridCellCohorts.size(); ii++)
            {
                sum += GridCellCohorts[ii].size();
            }

            return sum;
        }

    };
//
/** \brief IEnumerator for the grid cell cohorts */
//     class GridCellCohortsEnum : IEnumerator<List<Cohort>>
//    {
/** \brief The grid cell cohorts as a vector (with elements corresponding to functional groups) of lists of cohorts */
//         List<Cohort>[] GridCellCohorts;
//
/** \brief Current position in the vector of lists of cohorts */
//        int position = -1;
//
/** \brief Assign the passed set of grid cell cohorts to the internal vector of lists of cohorts
@param list  */
//         GridCellCohortsEnum(List<Cohort>[] list)
//        {
//            GridCellCohorts = list;
//        }
//
/** \brief Move to the next element in the vector of lists of cohorts
@return True if the end of the list had not been reached*/
//         bool MoveNext()
//        {
//            position++;
//            return (position <  GridCellCohorts.Length);
//        }
//
/** \brief Move back to the first element in the vector of lists of cohorts */
//         void Reset()
//        {
//            position = -1;
//        }
//
/** \brief Returns the list of cohorts for the current position (i.e. functional group) in the vector of lists of cohorts */
//        object System.Collections.IEnumerator.Current
//        {
//            get
//            {
//                return Current;
//            }
//        }
//
/** \brief Get the list of cohorts for the current position (i.e. functional group) in the vector of lists of cohorts */
//         List<Cohort> Current
//        {
//            get
//            {
//                try
//                {
//                    return GridCellCohorts[position];
//                }
//                catch (IndexOutOfRangeException)
//                {
//                    throw new InvalidOperationException();
//                }
//            }
//        }
//
/** \brief Destructor for the grid cell cohorts enumerator */
//         void Dispose()
//        {
//        }
//
//    }
//
//}
#endif
