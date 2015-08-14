#ifndef GRIDCELLSTOCKHANDLER_H
#define GRIDCELLSTOCKHANDLER_H
#include <Stock.h>
/** \file GridCellStockHandler.h
 * \brief the GridCellStockHandler header file
 */


//
//
//namespace Madingley
//{
/** \brief Handles the stocks in a grid cell */
//    /// <todoD>CAN USE COLLECTIONBASE SYNTAX TO ELIMINATE A LOT OF CODE HERE</todoD>
//    /// <todo>Create a wrapper class to handle our array of lists of gridCellStocks within an individual grid cell</todo>
class GridCellStockHandler //: IList<List<Stock>>, IEnumerable<List<Stock>>
    {
/** \brief A vector (with elements correpsonding to functional groups) of lists of stocks in the current grid cell*/
        vector< vector<Stock> > GridCellStocks;

    public:
vector<Stock>& operator[](unsigned i){return GridCellStocks[i];}
unsigned size(){return GridCellStocks.size();}
/** \brief Overloaded constructor for the grid cell stock handler: initialises a new vector of lists of stocks*/
        GridCellStockHandler()
       {
//            GridCellStocks = new List<Stock>[0];
       }

/** \brief Overloaded constructor for the grid cell stock handler: initialises a new vector of lists of stocks with number of elements equal to the number of functional groups 
@param NumFunctionalGroups The number of stock functional groups in the model */
//        GridCellStockHandler(int NumFunctionalGroups)
       void setSize(int NumFunctionalGroups)
       {
//            GridCellStocks = new List<Stock>[NumFunctionalGroups];
       }

/** \brief Overloaded constructor for the grid cell stock handler: update the grid cell stocks with the a set of existing stocks
@param ExistingStocks */  
//         GridCellStockHandler(List<Stock>[] ExistingStocks)
//        {
//            GridCellStocks = ExistingStocks;
//        }
//
/** \brief Get or set the list of stocks for a specified functional group index
@param index The functional group index 
@return The list of stocks from the specified functional group index */
//         List<Stock> this[int index]
//        {
//            get { return GridCellStocks[index]; }
//            set
//            {
//                GridCellStocks[index] = value;
//            }
//        }
//
//        // Gets of sets a Stock within the array of lists of gridCellStocks where the first element of the 2-element vector passed in is the array index and the second element is the list index
/** \brief Get or set the stock at a specified position within a specified functional group index
@param index Pair of values corresponding to the functional group index and the position of the stock within this functional group 
@return The stock at the specified position */
//         Stock this[int[] index]
//        {
//            get { return GridCellStocks[index[0]][index[1]]; }
//            set { GridCellStocks[index[0]][index[1]] = value; }
//        }
//
/** \brief Get the functional group index of the passed list of stocks
@param item The list of stocks to get the functional group index for 
@return The functional group index of the passed list of stocks */
//         int IndexOf(List<Stock> item)
//        {
//            return ((IList<List<Stock>>)GridCellStocks).IndexOf(item);
//        }
//
/** \brief NOT CURRENTLY USED
@param index NOT CURRENTLY USED 
@param item NOT CURRENTLY USED */
//         void Insert(int index, List<Stock> item)
//        {
//            throw new NotImplementedException();
//        }
//        
/** \brief NOT CURRENTLY USED
@param index NOT CURRENTLY USED */
//         void RemoveAt(int index)
//        {
//            throw new NotImplementedException();
//        }
//
/** \brief NOT CURRENTLY USED
@param item NOT CURRENTLY USED */
//         void Add(List<Stock> item)
//        {
//            throw new NotImplementedException();
//        }
//
/** \brief NOT CURRENTLY USED*/
//         void Clear()
//        {
//            throw new NotImplementedException();
//        }
//
/** \brief NOT CURRENTLY USED
@param item NOT CURRENTLY USED 
@return NOT CURRENTLY USED */
//         bool Contains(List<Stock> item)
//        {
//            throw new NotImplementedException();
//        }
//
/** \brief NOT CURRENTLY USED
@param array NOT CURRENTLY USED 
@param arrayIndex NOT CURRENTLY USED */
//         void CopyTo(List<Stock>[] array, int arrayIndex)
//        {
//            throw new NotImplementedException();
//        }
//
/** \brief Get the number of stock functional groups*/
//         int Count
//        {
//            get { return GridCellStocks.Count(); }
//        }
//
/** \brief NOT CURRENTLY USED*/
//         bool IsReadOnly
//        {
//            get { throw new NotImplementedException(); }
//        }
//
/** \brief NOT CURRENTLY USED
@param item NOT CURRENTLY USED 
@return NOT CURRENTLY USED*/
//         bool Remove(List<Stock> item)
//        {
//            throw new NotImplementedException();
//        }
//
/** \brief Return the grid cell stocks as an IEnumerator
@return */
//         IEnumerator<List<Stock>> GetEnumerator()
//        {
//            return new GridCellStocksEnum(GridCellStocks);
//        }
//
/** \brief Return an IEnumerable as an IEnumerator
@return */
//        IEnumerator IEnumerable.GetEnumerator()
//        {
//            return (IEnumerator)GetEnumerator();
//        }
    };
//
/** \brief //    /// IEnumerator for the grid cell stocks
//    /// </summary>
//     class GridCellStocksEnum : IEnumerator<List<Stock>>
//    {
//        // The array of lists of gridCellStocks
/** \brief The grid cell stocks as a vector (with elements corresponding to functional groups) of lists of stocks*/
//         List<Stock>[] GridCellStocks;
//
/** \brief Current position in the vector of lists of stocks*/
//        int position = -1;
//
/** \brief Assign the passed set of grid cell stocks to the internal vector of lists of stocks 
@param list  */
//         GridCellStocksEnum(List<Stock>[] list)
//        {
//            GridCellStocks = list;
//        }
//
/** \brief Move to the next element in the vector of lists of stocks
@return True if the end of the list had not been reached*/
//         bool MoveNext()
//        {
//            position++;
//            return (position < GridCellStocks.Length);
//        }
//
/** \brief Move back to the first element in the vector of lists of stocks*/
//         void Reset()
//        {
//            position = -1;
//        }
//
/** \brief Returns the list of stocks for the current position (i.e. functional group) in the vector of lists of stocks*/
//        object System.Collections.IEnumerator.Current
//        {
//            get
//            {
//                return Current;
//            }
//        }
//
/** \brief Get the list of stocks for the current position (i.e. functional group) in the vector of lists of stocks*/
//         List<Stock> Current
//        {
//            get
//            {
//                try
//                {
//                    return GridCellStocks[position];
//                }
//                catch (IndexOutOfRangeException)
//                {
//                    throw new InvalidOperationException();
//                }
//            }
//        }
//
/** \brief Destructor for the grid cell stocks enumerator*/
//         void Dispose()
//        {
//        }
//    }
//
//}
#endif
