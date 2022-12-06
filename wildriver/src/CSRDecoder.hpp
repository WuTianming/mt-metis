/**
* @file CSRDecoder.hpp
* @brief A class for reading in a matrix into a CSR structure.
* @author Dominique LaSalle <dominique@solidlake.com>
* Copyright 2017
* @version 1
* @date 2017-07-25
*/




#ifndef WILDRIVER_CSRDECODER_HPP
#define WILDRIVER_CSRDECODER_HPP




#include "IRowMatrixReader.hpp"




namespace WildRiver
{


class CSRDecoder
{
  public:
    /**
    * @brief Create a new CSRDecoder.
    *
    * @param reader The reader to use.
    */
    CSRDecoder(
      IRowMatrixReader * reader);


    /**
     * @brief Get the number of rows, columns, and non-zeros in the matrix.
     *
     * @param nrows The number of rows.
     * @param ncols The number of columns.
     * @param nnz THe number of non-zeros.
     */
    void getInfo(
        dim_t & nrows,
        dim_t & ncols,
        ind_t & nnz);


    /**
     * @brief Get the sparse matrix in CSR form. The pointers must be
     * pre-allocated to the sizes required by the info of the matrix
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
    void read(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval,
        double * progress);


  private:
    /**
    * @brief Whether or not the matrix information has been set.
    */
    bool m_infoSet;


    /**
     * @brief The number of rows in the matrix.
     */
    dim_t m_numRows;


    /**
     * @brief The number of columns in the matrix.
     */
    dim_t m_numCols;


    /**
     * @brief The number of non-zeros in the matrix.
     */
    ind_t m_nnz;


    /**
    * @brief The reader to get the matrix row by row.
    */
    IRowMatrixReader * m_reader;




};




}




#endif
