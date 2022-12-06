/**
* @file GraphMatrixReader.hpp
* @brief An adapter class for treating graphs as matrices.
* @author Dominique LaSalle <dominique@solidlake.com>
* Copyright 2017, Solid Lake LLC
* @version 1
* @date 2017-07-24
*/




#ifndef WILDRIVER_GRAPHMATRIXREADER_HPP
#define WILDRIVER_GRAPHMATRIXREADER_HPP




#include <vector>
#include <string>
#include <memory>


#include "IMatrixReader.hpp"
#include "IGraphReader.hpp"




namespace WildRiver
{


class GraphMatrixReader :
  public IMatrixReader
{
  public:
    /**
    * @brief Create a new adapter.
    *
    * @param reader The MatrixReader to adapt.
    */
    GraphMatrixReader(
        std::unique_ptr<IGraphReader>& reader);


    /**
    * @brief Destructor.
    */
    virtual ~GraphMatrixReader();


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
        double * progress) override;


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
        ind_t & nnz) override;


  private:
    std::unique_ptr<IGraphReader> m_reader;


    // disable copying
    GraphMatrixReader(
        GraphMatrixReader const & rhs);
    GraphMatrixReader& operator=(
        GraphMatrixReader const & rhs);




};




}




#endif
