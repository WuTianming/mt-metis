/**
 * @file IMatrixWriter.hpp
 * @brief Interface for writing out sparse matrices.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2015-09-02
 */




#ifndef WILDRIVER_IMATRIXWRITER_HPP
#define WILDRIVER_IMATRIXWRITER_HPP




#include <string>
#include <vector>
#include <stdexcept>

#include "base.h"




namespace WildRiver 
{


class IMatrixWriter
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IMatrixWriter() 
    {
      // do nothing
    }


    /**
     * @brief Set the matrix information for this file.
     *
     * @param nrows The number of rows in the matrix.
     * @param ncols The number of columns in the matrix.
     * @param nnz The number of non-zeroes in the matrix.
     */
    virtual void setInfo(
        dim_t nrows,
        dim_t ncols,
        ind_t nnz) = 0;


    /**
     * @brief Write the given CSR structure to the file. The information for
     * the matrix must already be set.
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     */
    virtual void write(
        ind_t const * rowptr,
        dim_t const * rowind,
        val_t const * rowval) = 0;


};


}




#endif
