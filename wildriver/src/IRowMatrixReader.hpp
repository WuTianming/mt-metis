/**
* @file IRowMatrixReader.hpp
* @brief An interface for reading in matrices row by row.
* @author Dominique LaSalle <dominique@solidlake.com>
* Copyright 2017, Solid Lake LLC
* @version 1
* @date 2017-07-25
*/




#ifndef WILDRIVER_IROWMATRIXREADER_HPP
#define WILDRIVER_IROWMATRIXREADER_HPP




#include "base.h"




namespace WildRiver
{


class IRowMatrixReader
{
  public:
    /**
    * @brief Destructor.
    */
    virtual ~IRowMatrixReader()
    {
      // do nothing
    }


    /**
     * @brief Read the header of this matrix file. Populates internal fields
     * with the header information.
     */
    virtual void readHeader(
        dim_t & numRows,
        dim_t & numCols,
        ind_t & nnz) = 0;


    /**
     * @brief Get the next row in the matrix (adjacecny list in the graph).
     *
     * @param numNonZeros The number of non-zeros in the row (output).
     * @param columns The column of each non-zero entry (must be of length at
     * least the number of non-zero entries).
     * @param values The value of each non-zero entry (must be null or of
     * length at least the number of non-zero entries).
     *
     * @return True if another row was found in the file.
     */
    virtual void getNextRow(
        dim_t * numNonZeros,
        dim_t * columns,
        val_t * values) = 0;




};




}




#endif
