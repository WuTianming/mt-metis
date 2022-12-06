/**
 * @file IMatrixReader.hpp
 * @brief Interface for reading sparse matrices.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_IMATRIXREADER_HPP
#define WILDRIVER_IMATRIXREADER_HPP




#include <string>
#include <vector>

#include "base.h"




namespace WildRiver 
{


class IMatrixReader
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IMatrixReader() 
    {
      // do nothing
    }


    /**
     * @brief Get the matrix in CSR form. The pointers must be
     * pre-allocated to the sizes required by the info of the matrix.
     *
     * |rowptr| = nrows + 1
     * |rowind| = nnz
     * |rowval| = nnz
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     * @param progress The variable to update as the matrix is loaded (may be
     * null).
     */
    virtual void read(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval,
        double * progress) = 0;


    /**
     * @brief Get the number of rows, columns, and non-zeros in the matrix.
     *
     * @param nrows The number of rows.
     * @param ncols The number of columns.
     * @param nnz THe number of non-zeros.
     */
    virtual void getInfo(
        dim_t & nrows,
        dim_t & ncols,
        ind_t & nnz) = 0;


};




}




#endif
