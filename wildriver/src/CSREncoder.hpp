/**
 * @file CSREncoder.hpp
 * @brief The CSREncoder class.
 * @author Dominique LaSalle <dominique@solidlake.com>
 * Copyright 2017-2018
 * @version 1
 * @date 2017-07-26
 */




#ifndef WILDRIVER_CSRENCODER_HPP
#define WILDRIVER_CSRENCODER_HPP




#include "IRowMatrixWriter.hpp"
#include <memory>




namespace WildRiver
{


/**
* @brief A class for writing matrix row by row from a CSR data structure
*/
class CSREncoder
{
  public:
    /**
     * @brief Create a new CSREncoder.
     *
     * @param writer The writer to use.
     */
    CSREncoder(
        IRowMatrixWriter * writer);


    /**
     * @brief Set the matrix information.
     *
     * @param nrows The number of rows in the matrix.
     * @param ncols The number of columns in the matrix.
     * @param nnz The number of non-zeroes in the matrix.
     */
    void setInfo(
        dim_t nrows,
        dim_t ncols,
        ind_t nnz);


    /**
     * @brief Write the given CSR structure. The information for
     * the matrix must already be set.
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     */
    void write(
        ind_t const * rowptr,
        dim_t const * rowind,
        val_t const * rowval);




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
     * @brief The writer to set the matrix row by row.
     */
    IRowMatrixWriter * m_writer;




};




}




#endif
